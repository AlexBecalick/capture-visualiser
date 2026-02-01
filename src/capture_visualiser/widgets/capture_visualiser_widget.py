"""Main Capture Visualiser widget."""

from __future__ import annotations

from qtpy.QtCore import Qt, QThread, QTimer
from qtpy.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QProgressBar,
    QPushButton,
    QVBoxLayout,
    QWidget,
)
from napari.viewer import Viewer

from capture_visualiser._atlases import get_available_atlases
from capture_visualiser._constants import (
    DEFAULT_ATLAS_NAME,
    DEFAULT_CAPTURE_HEIGHT_MM,
    DEFAULT_CAPTURE_WIDTH_MM,
    XENIUM_SLIDE_HEIGHT_MM,
    XENIUM_SLIDE_WIDTH_MM,
)
from capture_visualiser._presets import get_presets
from capture_visualiser.atlas_layer import PLANE_CONFIG, AtlasCoronalView
from capture_visualiser.atlas_loader import AtlasInstallWorker, is_atlas_downloaded
from capture_visualiser.capture_rectangle import (
    CAPTURE_COLORS,
    SHAPES_PER_CAPTURE,
    add_capture_shape_to_layer_3d,
    add_default_capture_for_plane,
    get_shapes_physical_params,
    infer_shape_face_plane,
    recreate_shapes_in_new_atlas,
    rotate_capture_shape_at_index_3d,
    update_capture_shape_at_index_3d,
    update_capture_slice_positions,
)
from capture_visualiser.widgets.ap_slider import APSliderWidget


class CaptureVisualiserWidget(QWidget):
    """Main widget for Capture Visualiser: atlas + AP slider + capture rectangle."""

    def __init__(self, napari_viewer: Viewer, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._viewer = napari_viewer
        self._atlas_view: AtlasCoronalView | None = None
        self._ap_slider: APSliderWidget | None = None
        self._ref_layer_name: str | None = None
        self._ann_layer_name: str | None = None
        self._capture_layer_name: str | None = None
        self._capture_color_index: int = 0
        self._shape_colors: list[str] = []
        self._install_thread: QThread | None = None
        self._install_worker: AtlasInstallWorker | None = None
        self._canonical_params: list[dict] | None = None
        self._selection_sync_block = False
        self._last_capture_slice: int | None = None
        self._updating_capture_positions = False

        layout = QVBoxLayout(self)

        # Header
        header = QLabel(
            "<b>Capture Visualiser</b><br>"
            "Overlay spatial transcriptomics capture areas onto brain atlas coronal sections."
        )
        header.setWordWrap(True)
        layout.addWidget(header)

        # Atlas section
        atlas_group = QGroupBox("Atlas")
        atlas_layout = QVBoxLayout(atlas_group)
        atlas_row = QHBoxLayout()
        self._atlas_combo = QComboBox()
        self._atlas_combo.setToolTip("Select a BrainGlobe atlas. Downloads on first use if needed.")
        for atlas_id, display_name in get_available_atlases():
            self._atlas_combo.addItem(display_name, atlas_id)
        # Set default to allen_mouse_25um
        default_idx = self._atlas_combo.findData(DEFAULT_ATLAS_NAME)
        if default_idx >= 0:
            self._atlas_combo.setCurrentIndex(default_idx)
        atlas_row.addWidget(self._atlas_combo)
        self._load_btn = QPushButton("Load atlas")
        self._load_btn.setToolTip("Load the selected atlas. Downloads on first use if needed.")
        self._load_btn.clicked.connect(self._on_load_atlas)
        atlas_row.addWidget(self._load_btn)
        atlas_layout.addLayout(atlas_row)
        self._progress_bar = QProgressBar()
        self._progress_bar.setVisible(False)
        self._progress_bar.setRange(0, 100)
        atlas_layout.addWidget(self._progress_bar)
        layout.addWidget(atlas_group)

        # Preset buttons
        presets_group = QGroupBox("Presets")
        presets_layout = QVBoxLayout(presets_group)
        for preset in get_presets():
            name = preset["name"]
            desc = preset["description"]
            w, h = preset["width_mm"], preset["height_mm"]
            btn = QPushButton(f"{name} ({w} x {h} mm)")
            btn.setToolTip(desc)
            btn.clicked.connect(
                lambda checked=False, pw=w, ph=h: self._apply_preset(pw, ph)
            )
            presets_layout.addWidget(btn)
        layout.addWidget(presets_group)

        # Custom capture area size
        capture_group = QGroupBox("Custom capture area")
        capture_group.setToolTip(
            "Resize the selected capture area. Select a shape first with the Select tool."
        )
        capture_layout = QFormLayout(capture_group)
        self._width_spin = QDoubleSpinBox()
        self._width_spin.setRange(0.5, 100.0)
        self._width_spin.setValue(DEFAULT_CAPTURE_WIDTH_MM)
        self._width_spin.setSuffix(" mm")
        self._width_spin.setDecimals(1)
        capture_layout.addRow("Width:", self._width_spin)

        self._height_spin = QDoubleSpinBox()
        self._height_spin.setRange(0.5, 100.0)
        self._height_spin.setValue(DEFAULT_CAPTURE_HEIGHT_MM)
        self._height_spin.setSuffix(" mm")
        self._height_spin.setDecimals(1)
        capture_layout.addRow("Height:", self._height_spin)
        self._apply_custom_btn = QPushButton("Apply")
        self._apply_custom_btn.setToolTip("Resize the selected capture area to the dimensions above.")
        self._apply_custom_btn.clicked.connect(self._on_apply_custom)
        capture_layout.addRow("", self._apply_custom_btn)
        layout.addWidget(capture_group)

        self._rotate_btn = QPushButton("Rotate current capture area 90°")
        self._rotate_btn.setToolTip("Rotate the selected capture area by 90 degrees clockwise.")
        self._rotate_btn.clicked.connect(self._on_rotate_90)
        self._rotate_btn.setEnabled(False)
        layout.addWidget(self._rotate_btn)

        self._create_new_btn = QPushButton("Create a new capture area")
        self._create_new_btn.setToolTip(
            "Add another capture area. Uses the dimensions above. Select shapes to resize them."
        )
        self._create_new_btn.clicked.connect(self._on_create_new_capture_area)
        self._create_new_btn.setEnabled(False)
        layout.addWidget(self._create_new_btn)

        # Plane and slice slider (created when atlas loads)
        self._ap_group = QGroupBox("View plane")
        self._ap_placeholder = QLabel("Load atlas first to browse planes.")
        ap_layout = QVBoxLayout(self._ap_group)
        ap_layout.addWidget(self._ap_placeholder)
        plane_row = QHBoxLayout()
        self._plane_coronal_btn = QPushButton("Coronal")
        self._plane_coronal_btn.setCheckable(True)
        self._plane_coronal_btn.setToolTip("View coronal sections (AP slice, DV × ML displayed)")
        self._plane_coronal_btn.clicked.connect(lambda: self._on_plane_changed("coronal"))
        self._plane_sagittal_btn = QPushButton("Sagittal")
        self._plane_sagittal_btn.setCheckable(True)
        self._plane_sagittal_btn.setToolTip("View sagittal sections (ML slice, AP × DV displayed)")
        self._plane_sagittal_btn.clicked.connect(lambda: self._on_plane_changed("sagittal"))
        self._plane_horizontal_btn = QPushButton("Horizontal")
        self._plane_horizontal_btn.setCheckable(True)
        self._plane_horizontal_btn.setToolTip("View horizontal sections (DV slice, AP × ML displayed)")
        self._plane_horizontal_btn.clicked.connect(lambda: self._on_plane_changed("horizontal"))
        for btn in (self._plane_coronal_btn, self._plane_sagittal_btn, self._plane_horizontal_btn):
            btn.setEnabled(False)
            plane_row.addWidget(btn)
        ap_layout.addLayout(plane_row)
        layout.addWidget(self._ap_group)

        # Instructions
        instructions = QLabel(
            "<i>Use napari's Select shapes tool to select, move and rotate capture areas. "
            "Click a preset or Apply to resize the selected area.</i>"
        )
        instructions.setWordWrap(True)
        layout.addWidget(instructions)

        layout.addStretch()

    def _on_load_atlas(self) -> None:
        atlas_name = self._atlas_combo.currentData()
        if not atlas_name:
            atlas_name = DEFAULT_ATLAS_NAME

        if not is_atlas_downloaded(atlas_name):
            self._download_then_load(atlas_name)
            return

        self._do_load_atlas(atlas_name)

    def _download_then_load(self, atlas_name: str) -> None:
        """Download atlas with progress bar, then load."""
        self._load_btn.setEnabled(False)
        self._atlas_combo.setEnabled(False)
        self._progress_bar.setVisible(True)
        self._progress_bar.setRange(0, 0)  # Indeterminate until we get total

        self._install_thread = QThread()
        self._install_worker = AtlasInstallWorker(atlas_name)
        self._install_worker.moveToThread(self._install_thread)

        def on_progress(completed: int, total: int) -> None:
            if total > 0:
                self._progress_bar.setRange(0, total)
                self._progress_bar.setValue(completed)

        def on_finished() -> None:
            self._install_thread.quit()
            self._progress_bar.setVisible(False)
            self._load_btn.setEnabled(True)
            self._atlas_combo.setEnabled(True)
            self._do_load_atlas(atlas_name)

        def on_error(msg: str) -> None:
            self._install_thread.quit()
            self._progress_bar.setVisible(False)
            self._load_btn.setEnabled(True)
            self._atlas_combo.setEnabled(True)
            self._show_error(f"Atlas download failed: {msg}")

        self._install_worker.progress.connect(on_progress)
        self._install_worker.finished.connect(on_finished)
        self._install_worker.error.connect(on_error)
        self._install_thread.started.connect(self._install_worker.run)
        self._install_thread.start()

    def _do_load_atlas(self, atlas_name: str) -> None:
        """Load atlas and add to viewer. Replaces existing atlas if any."""
        # Save existing capture shapes for atlas switch
        saved_params = None
        saved_plane = "coronal"
        if (
            self._atlas_view is not None
            and self._capture_layer_name is not None
            and self._capture_layer_name in self._viewer.layers
        ):
            saved_plane = self._atlas_view.plane
            layer = self._viewer.layers[self._capture_layer_name]
            res_per_axis = (
                self._atlas_view.res_ap,
                self._atlas_view.res_dv,
                self._atlas_view.res_ml,
            )
            slice_axis = PLANE_CONFIG[self._atlas_view.plane]["slice_axis"]
            slice_pos_um = self._atlas_view.current_slice_index() * res_per_axis[slice_axis]
            saved_params = get_shapes_physical_params(
                layer, self._atlas_view.plane,
                self._atlas_view.res_ap, self._atlas_view.res_dv, self._atlas_view.res_ml,
                slice_position_um=slice_pos_um,
                prev_params=self._canonical_params,
            )

        # Remove old atlas layers
        if self._ref_layer_name and self._ref_layer_name in self._viewer.layers:
            self._viewer.layers.remove(self._ref_layer_name)
        if self._ann_layer_name and self._ann_layer_name in self._viewer.layers:
            self._viewer.layers.remove(self._ann_layer_name)

        try:
            self._atlas_view = AtlasCoronalView(atlas_name)
        except Exception as e:
            self._show_error(f"Failed to load atlas: {e}")
            return

        self._ref_layer_name, self._ann_layer_name = self._atlas_view.add_to_viewer(
            self._viewer
        )

        if saved_params is not None and self._capture_layer_name and self._capture_layer_name in self._viewer.layers:
            self._canonical_params = saved_params
            layer = self._viewer.layers[self._capture_layer_name]
            recreate_shapes_in_new_atlas(
                layer,
                saved_params,
                saved_plane,
                self._atlas_view.res_ap,
                self._atlas_view.res_dv,
                self._atlas_view.res_ml,
                self._shape_colors,
            )
            self._atlas_view.apply_plane(self._viewer, saved_plane)
        else:
            self._canonical_params = None
            # Initial load: add first capture rectangle
            self._capture_color_index = 0
            self._shape_colors = [CAPTURE_COLORS[0]]
            self._capture_layer_name = add_default_capture_for_plane(
                self._viewer,
                "coronal",
                self._atlas_view.res_ap,
                self._atlas_view.res_dv,
                self._atlas_view.res_ml,
                self._atlas_view.n_ap_slices,
                self._atlas_view.n_dv_slices,
                self._atlas_view.n_ml_slices,
                slice_index=self._atlas_view.ap_index,
            )
            self._width_spin.setValue(XENIUM_SLIDE_WIDTH_MM)
            self._height_spin.setValue(XENIUM_SLIDE_HEIGHT_MM)
            layer = self._viewer.layers[self._capture_layer_name]
            layer.selected_data = {0}
            slice_pos_um = self._atlas_view.ap_index * self._atlas_view.res_ap
            self._canonical_params = get_shapes_physical_params(
                layer, "coronal",
                self._atlas_view.res_ap, self._atlas_view.res_dv, self._atlas_view.res_ml,
                slice_position_um=slice_pos_um,
            )

        # Move capture layer above both atlas layers (reference and annotation)
        cap_layer = self._viewer.layers[self._capture_layer_name]
        cap_idx = self._viewer.layers.index(cap_layer)
        self._viewer.layers.move(cap_idx, len(self._viewer.layers))
        self._connect_capture_selection_sync()
        # Re-assert dims range (2D shapes layer can corrupt napari's computed range)
        self._atlas_view.set_dims_range(self._viewer)
        self._viewer.dims.current_step = tuple(self._atlas_view._slice_indices)
        QTimer.singleShot(0, self._reassert_dims_range)

        # AP slider
        if self._ap_placeholder is not None:
            self._ap_group.layout().removeWidget(self._ap_placeholder)
            self._ap_placeholder.deleteLater()
            self._ap_placeholder = None

        if self._ap_slider is not None:
            self._ap_group.layout().removeWidget(self._ap_slider)
        plane = saved_plane if saved_params else "coronal"
        cfg = PLANE_CONFIG[plane]
        label = f"{cfg['slice_label']} plane:"
        self._ap_slider = APSliderWidget(
            min_val=0,
            max_val=self._atlas_view.n_slices_for_plane() - 1,
            value=self._atlas_view.current_slice_index(),
        )
        self._update_plane_buttons_state()
        self._ap_slider.set_label(label)
        self._ap_slider.valueChanged.connect(self._on_slice_slider_changed)
        self._ap_group.layout().addWidget(self._ap_slider)
        self._update_plane_buttons_state()
        for btn in (self._plane_coronal_btn, self._plane_sagittal_btn, self._plane_horizontal_btn):
            btn.setEnabled(True)
        self._update_create_rotate_for_plane()

        try:
            self._viewer.dims.events.current_step.disconnect(self._on_napari_dims_changed)
        except (TypeError, ValueError):
            pass
        self._viewer.dims.events.current_step.connect(self._on_napari_dims_changed)

        self._create_new_btn.setEnabled(True)
        self._rotate_btn.setEnabled(True)

    def _update_plane_buttons_state(self) -> None:
        """Update check state of plane buttons to reflect current plane."""
        plane = self._atlas_view.plane if self._atlas_view else "coronal"
        self._plane_coronal_btn.setChecked(plane == "coronal")
        self._plane_sagittal_btn.setChecked(plane == "sagittal")
        self._plane_horizontal_btn.setChecked(plane == "horizontal")

    def _reassert_dims_range(self) -> None:
        """Re-assert atlas dims range and center (called deferred to override async napari updates)."""
        if self._atlas_view is not None:
            self._atlas_view.reset_to_center()
            self._atlas_view.set_dims_range(self._viewer)
            self._viewer.dims.current_step = tuple(self._atlas_view._slice_indices)

    def _update_create_rotate_for_plane(self) -> None:
        """Enable create/rotate when atlas is loaded."""
        self._create_new_btn.setEnabled(self._atlas_view is not None)
        self._rotate_btn.setEnabled(self._atlas_view is not None)

    def _on_plane_changed(self, plane: str) -> None:
        """Switch view to the selected plane. Keeps existing capture areas (they retain their orientation)."""
        if self._atlas_view is None:
            return
        if plane == self._atlas_view.plane:
            return
        self._last_capture_slice = None
        self._atlas_view.reset_to_center()
        self._atlas_view.apply_plane(self._viewer, plane)
        if self._capture_layer_name and self._capture_layer_name in self._viewer.layers:
            cap_layer = self._viewer.layers[self._capture_layer_name]
            cap_idx = self._viewer.layers.index(cap_layer)
            self._viewer.layers.move(cap_idx, len(self._viewer.layers))
        self._atlas_view.set_dims_range(self._viewer)
        self._viewer.dims.current_step = tuple(self._atlas_view._slice_indices)
        QTimer.singleShot(0, self._reassert_dims_range)
        cfg = PLANE_CONFIG[plane]
        self._ap_slider.set_label(f"{cfg['slice_label']} plane:")
        self._ap_slider.set_range(0, self._atlas_view.n_slices_for_plane() - 1)
        self._ap_slider.set_value(self._atlas_view.current_slice_index())
        self._update_plane_buttons_state()
        self._update_create_rotate_for_plane()

    def _on_slice_slider_changed(self, value: int) -> None:
        if self._atlas_view is None:
            return
        self._atlas_view.set_current_slice_index(value)
        step = tuple(self._atlas_view._slice_indices)
        self._viewer.dims.current_step = step

    def _get_selected_shape_index(self) -> int | None:
        """Return a selected shape index (first if multiple), or None if none selected."""
        if (
            self._capture_layer_name is None
            or self._capture_layer_name not in self._viewer.layers
        ):
            return None
        layer = self._viewer.layers[self._capture_layer_name]
        sel = layer.selected_data
        if len(sel) < 1:
            return None
        return next(iter(sel))

    def _expand_selection_to_cuboid(self) -> None:
        """When user selects a shape, expand selection to include all faces of that cuboid."""
        if (
            self._selection_sync_block
            or not self._capture_layer_name
            or self._capture_layer_name not in self._viewer.layers
        ):
            return
        layer = self._viewer.layers[self._capture_layer_name]
        sel = set(layer.selected_data)
        if not sel:
            return
        expanded = set()
        for idx in sel:
            group = self._get_capture_group_indices(idx)
            expanded.update(group)
        if expanded != sel:
            self._selection_sync_block = True
            try:
                layer.selected_data = expanded
            finally:
                self._selection_sync_block = False

    def _connect_capture_selection_sync(self) -> None:
        """Connect selection-changed callback so cuboid faces move together."""
        if not self._capture_layer_name or self._capture_layer_name not in self._viewer.layers:
            return
        layer = self._viewer.layers[self._capture_layer_name]
        try:
            layer.selected_data.events.items_changed.disconnect(self._expand_selection_to_cuboid)
        except (TypeError, ValueError):
            pass
        layer.selected_data.events.items_changed.connect(self._expand_selection_to_cuboid)

    def _get_capture_group_indices(self, shape_index: int) -> list[int]:
        """Return all shape indices for the capture area containing shape_index."""
        group_start = (shape_index // SHAPES_PER_CAPTURE) * SHAPES_PER_CAPTURE
        layer = self._viewer.layers[self._capture_layer_name]
        n = layer.nshapes
        return [
            group_start + j
            for j in range(SHAPES_PER_CAPTURE)
            if group_start + j < n
        ]

    def _apply_preset(self, width_mm: float, height_mm: float) -> None:
        """Apply a preset: resize the selected capture area and update spinboxes."""
        idx = self._get_selected_shape_index()
        if idx is None:
            self._show_error("Select a capture area first (Select shapes tool).")
            return
        self._width_spin.blockSignals(True)
        self._height_spin.blockSignals(True)
        self._width_spin.setValue(width_mm)
        self._height_spin.setValue(height_mm)
        self._width_spin.blockSignals(False)
        self._height_spin.blockSignals(False)
        self._resize_shape_at_index(idx, width_mm, height_mm)

    def _on_apply_custom(self) -> None:
        """Apply custom dimensions to the selected capture area."""
        idx = self._get_selected_shape_index()
        if idx is None:
            self._show_error("Select a capture area first (Select shapes tool).")
            return
        w = self._width_spin.value()
        h = self._height_spin.value()
        self._resize_shape_at_index(idx, w, h)

    def _resize_shape_at_index(
        self, shape_index: int, width_mm: float, height_mm: float
    ) -> None:
        if self._atlas_view is None or self._capture_layer_name is None:
            return
        if self._capture_layer_name not in self._viewer.layers:
            return
        layer = self._viewer.layers[self._capture_layer_name]
        shape_indices = self._get_capture_group_indices(shape_index)
        plane = infer_shape_face_plane(layer.data[shape_index])
        update_capture_shape_at_index_3d(
            layer,
            shape_indices=shape_indices,
            width_mm=width_mm,
            height_mm=height_mm,
            plane=plane,
            res_ap_um=self._atlas_view.res_ap,
            res_dv_um=self._atlas_view.res_dv,
            res_ml_um=self._atlas_view.res_ml,
        )
        self._atlas_view.set_dims_range(self._viewer)
        self._viewer.dims.current_step = tuple(self._atlas_view._slice_indices)

    def _on_create_new_capture_area(self) -> None:
        """Add a new capture area with current dimensions and next color."""
        if (
            self._atlas_view is None
            or self._capture_layer_name is None
            or self._capture_layer_name not in self._viewer.layers
        ):
            self._show_error("Load the atlas first.")
            return
        layer = self._viewer.layers[self._capture_layer_name]
        self._capture_color_index = (self._capture_color_index + 1) % len(CAPTURE_COLORS)
        color = CAPTURE_COLORS[self._capture_color_index]
        self._shape_colors.append(color)
        slice_idx = self._atlas_view.current_slice_index()
        all_face_colors = [
            c for c in self._shape_colors
            for _ in range(SHAPES_PER_CAPTURE)
        ]
        add_capture_shape_to_layer_3d(
            layer,
            width_mm=self._width_spin.value(),
            height_mm=self._height_spin.value(),
            plane=self._atlas_view.plane,
            res_ap_um=self._atlas_view.res_ap,
            res_dv_um=self._atlas_view.res_dv,
            res_ml_um=self._atlas_view.res_ml,
            face_color=color,
            all_face_colors=all_face_colors,
            shape_ap=self._atlas_view.n_ap_slices,
            shape_dv=self._atlas_view.n_dv_slices,
            shape_ml=self._atlas_view.n_ml_slices,
            slice_index=slice_idx,
        )
        layer.selected_data = {layer.nshapes - 1}
        res = (self._atlas_view.res_ap, self._atlas_view.res_dv, self._atlas_view.res_ml)
        slice_axis = PLANE_CONFIG[self._atlas_view.plane]["slice_axis"]
        slice_pos_um = slice_idx * res[slice_axis]
        self._canonical_params = get_shapes_physical_params(
            layer, self._atlas_view.plane,
            self._atlas_view.res_ap, self._atlas_view.res_dv, self._atlas_view.res_ml,
            slice_position_um=slice_pos_um,
            prev_params=self._canonical_params,
        )

    def _on_rotate_90(self) -> None:
        """Rotate the selected capture area by 90 degrees."""
        idx = self._get_selected_shape_index()
        if idx is None:
            self._show_error("Select a capture area first (Select shapes tool).")
            return
        if self._atlas_view is None or self._capture_layer_name is None:
            return
        if self._capture_layer_name not in self._viewer.layers:
            return
        layer = self._viewer.layers[self._capture_layer_name]
        shape_indices = self._get_capture_group_indices(idx)
        plane = infer_shape_face_plane(layer.data[idx])
        rotate_capture_shape_at_index_3d(
            layer,
            shape_indices=shape_indices,
            angle_deg=90.0,
            plane=plane,
            res_ap_um=self._atlas_view.res_ap,
            res_dv_um=self._atlas_view.res_dv,
            res_ml_um=self._atlas_view.res_ml,
        )

    def _on_napari_dims_changed(self) -> None:
        if self._ap_slider is None or self._atlas_view is None:
            return
        step = self._viewer.dims.current_step
        self._atlas_view.sync_step_from_viewer(step)
        current_slice = self._atlas_view.current_slice_index()
        if self._ap_slider.value() != current_slice:
            self._ap_slider.blockSignals(True)
            self._ap_slider.set_value(current_slice)
            self._ap_slider.blockSignals(False)
        if (
            self._capture_layer_name
            and self._capture_layer_name in self._viewer.layers
            and self._last_capture_slice != current_slice
            and not self._updating_capture_positions
        ):
            self._last_capture_slice = current_slice
            plane = self._atlas_view.plane
            layer = self._viewer.layers[self._capture_layer_name]
            slice_val = float(current_slice)

            def _do_update(pl=plane, ly=layer, sv=slice_val) -> None:
                if self._updating_capture_positions:
                    return
                self._updating_capture_positions = True
                try:
                    update_capture_slice_positions(ly, pl, sv)
                finally:
                    self._updating_capture_positions = False

            QTimer.singleShot(0, _do_update)

    def _show_error(self, message: str) -> None:
        from napari.utils.notifications import show_error
        show_error(message)

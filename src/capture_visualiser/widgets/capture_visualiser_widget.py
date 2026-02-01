"""Main Capture Visualiser widget."""

from __future__ import annotations

from qtpy.QtWidgets import (
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)
from napari.viewer import Viewer

from capture_visualiser._constants import (
    DEFAULT_ATLAS_NAME,
    DEFAULT_CAPTURE_HEIGHT_MM,
    DEFAULT_CAPTURE_WIDTH_MM,
    XENIUM_SLIDE_HEIGHT_MM,
    XENIUM_SLIDE_WIDTH_MM,
)
from capture_visualiser._presets import get_presets
from capture_visualiser.atlas_layer import AtlasCoronalView
from capture_visualiser.capture_rectangle import (
    CAPTURE_COLORS,
    add_capture_rectangle,
    add_capture_shape_to_layer,
    rotate_capture_shape_at_index,
    update_capture_shape_at_index,
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
        self._load_btn = QPushButton("Load Allen Mouse 25 μm atlas")
        self._load_btn.setToolTip(
            "Load the Allen mouse brain atlas. Downloads on first use if needed."
        )
        self._load_btn.clicked.connect(self._on_load_atlas)
        atlas_layout.addWidget(self._load_btn)
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
        capture_layout.addRow("Width (ML):", self._width_spin)

        self._height_spin = QDoubleSpinBox()
        self._height_spin.setRange(0.5, 100.0)
        self._height_spin.setValue(DEFAULT_CAPTURE_HEIGHT_MM)
        self._height_spin.setSuffix(" mm")
        self._height_spin.setDecimals(1)
        capture_layout.addRow("Height (DV):", self._height_spin)
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

        # AP slider (created when atlas loads)
        self._ap_group = QGroupBox("AP plane")
        self._ap_placeholder = QLabel("Load atlas first to browse AP planes.")
        ap_layout = QVBoxLayout(self._ap_group)
        ap_layout.addWidget(self._ap_placeholder)
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
        try:
            atlas_name = DEFAULT_ATLAS_NAME
            self._atlas_view = AtlasCoronalView(atlas_name)
        except Exception as e:
            self._show_error(f"Failed to load atlas: {e}")
            return

        # Add to viewer
        self._ref_layer_name, self._ann_layer_name = self._atlas_view.add_to_viewer(
            self._viewer
        )

        # Add initial capture rectangle (Xenium default)
        self._capture_color_index = 0
        self._shape_colors = [CAPTURE_COLORS[0]]
        self._capture_layer_name = add_capture_rectangle(
            self._viewer,
            width_mm=XENIUM_SLIDE_WIDTH_MM,
            height_mm=XENIUM_SLIDE_HEIGHT_MM,
            res_dv_um=self._atlas_view.res_dv,
            res_ml_um=self._atlas_view.res_ml,
            shape_dv=self._atlas_view.shape_dv,
            shape_ml=self._atlas_view.shape_ml,
            face_color=self._shape_colors[0],
        )
        self._width_spin.setValue(XENIUM_SLIDE_WIDTH_MM)
        self._height_spin.setValue(XENIUM_SLIDE_HEIGHT_MM)
        # Select the initial shape so presets work immediately
        layer = self._viewer.layers[self._capture_layer_name]
        layer.selected_data = {0}

        # Replace AP placeholder with slider
        self._ap_group.layout().removeWidget(self._ap_placeholder)
        self._ap_placeholder.deleteLater()
        self._ap_placeholder = None

        self._ap_slider = APSliderWidget(
            min_val=0,
            max_val=self._atlas_view.n_ap_slices - 1,
            value=self._atlas_view.ap_index,
        )
        self._ap_slider.valueChanged.connect(self._on_ap_slider_changed)
        self._ap_group.layout().addWidget(self._ap_slider)

        # Sync napari dims when slider changes
        self._viewer.dims.events.current_step.connect(self._on_napari_dims_changed)

        # Disable load button, enable create new
        self._load_btn.setEnabled(False)
        self._load_btn.setText("Atlas loaded")
        self._create_new_btn.setEnabled(True)
        self._rotate_btn.setEnabled(True)

    def _on_ap_slider_changed(self, value: int) -> None:
        if self._atlas_view is None:
            return
        self._atlas_view.set_ap_index(value)
        step = list(self._viewer.dims.current_step)
        step[0] = value
        self._viewer.dims.current_step = tuple(step)

    def _get_selected_shape_index(self) -> int | None:
        """Return the single selected shape index, or None if not exactly one selected."""
        if (
            self._capture_layer_name is None
            or self._capture_layer_name not in self._viewer.layers
        ):
            return None
        layer = self._viewer.layers[self._capture_layer_name]
        sel = layer.selected_data
        if len(sel) != 1:
            return None
        return next(iter(sel))

    def _apply_preset(self, width_mm: float, height_mm: float) -> None:
        """Apply a preset: resize the selected capture area and update spinboxes."""
        idx = self._get_selected_shape_index()
        if idx is None:
            self._show_error("Select exactly one capture area first (Select shapes tool).")
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
            self._show_error("Select exactly one capture area first (Select shapes tool).")
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
        update_capture_shape_at_index(
            layer,
            shape_index=shape_index,
            width_mm=width_mm,
            height_mm=height_mm,
            res_dv_um=self._atlas_view.res_dv,
            res_ml_um=self._atlas_view.res_ml,
        )

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
        add_capture_shape_to_layer(
            layer,
            width_mm=self._width_spin.value(),
            height_mm=self._height_spin.value(),
            res_dv_um=self._atlas_view.res_dv,
            res_ml_um=self._atlas_view.res_ml,
            face_color=color,
            all_face_colors=self._shape_colors,
            shape_dv=self._atlas_view.shape_dv,
            shape_ml=self._atlas_view.shape_ml,
        )
        # Select the new shape
        layer.selected_data = {layer.nshapes - 1}

    def _on_rotate_90(self) -> None:
        """Rotate the selected capture area by 90 degrees."""
        idx = self._get_selected_shape_index()
        if idx is None:
            self._show_error("Select exactly one capture area first (Select shapes tool).")
            return
        if self._atlas_view is None or self._capture_layer_name is None:
            return
        if self._capture_layer_name not in self._viewer.layers:
            return
        layer = self._viewer.layers[self._capture_layer_name]
        rotate_capture_shape_at_index(
            layer,
            shape_index=idx,
            angle_deg=90.0,
            res_dv_um=self._atlas_view.res_dv,
            res_ml_um=self._atlas_view.res_ml,
        )

    def _on_napari_dims_changed(self) -> None:
        if self._ap_slider is None:
            return
        ap_step = int(self._viewer.dims.current_step[0])
        if self._ap_slider.value() != ap_step:
            self._ap_slider.blockSignals(True)
            self._ap_slider.set_value(ap_step)
            self._ap_slider.blockSignals(False)

    def _show_error(self, message: str) -> None:
        from napari.utils.notifications import show_error
        show_error(message)

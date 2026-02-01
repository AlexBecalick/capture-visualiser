"""Atlas view and layer management for coronal, sagittal, horizontal planes."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
from brainglobe_atlasapi import BrainGlobeAtlas

if TYPE_CHECKING:
    from napari.viewer import Viewer

PlaneMode = Literal["coronal", "sagittal", "horizontal"]

# dims.order: (slice_dim, disp_axis_0, disp_axis_1) -> slice first, display last two
# Atlas axes: 0=AP, 1=DV, 2=ML
PLANE_CONFIG = {
    "coronal": {
        "order": (0, 1, 2),   # slice AP, display (DV, ML)
        "slice_axis": 0,
        "slice_label": "AP",
        "display_axes": (1, 2),  # DV, ML
        "res": lambda r: (r[1], r[2]),  # res_dv, res_ml
    },
    "sagittal": {
        "order": (2, 0, 1),   # slice ML, display (AP, DV)
        "slice_axis": 2,
        "slice_label": "ML",
        "display_axes": (0, 1),
        "res": lambda r: (r[0], r[1]),  # res_ap, res_dv
    },
    "horizontal": {
        "order": (1, 0, 2),   # slice DV, display (AP, ML)
        "slice_axis": 1,
        "slice_label": "DV",
        "display_axes": (0, 2),
        "res": lambda r: (r[0], r[2]),  # res_ap, res_ml
    },
}


class AtlasCoronalView:
    """Manages brain atlas display for coronal, sagittal, horizontal planes.

    Atlas axes are (AP, DV, ML): axis 0 = AP, axis 1 = DV, axis 2 = ML.
    """

    def __init__(self, atlas_name: str) -> None:
        self.atlas = BrainGlobeAtlas(atlas_name=atlas_name)
        self._ref = np.asarray(self.atlas.reference)
        self._ann = np.asarray(self.atlas.annotation)
        self.resolution = self.atlas.resolution  # (res_ap, res_dv, res_ml) μm
        self._plane: PlaneMode = "coronal"
        # Slice indices for each axis
        self._slice_indices = [
            self._ref.shape[0] // 2,  # AP
            self._ref.shape[1] // 2,  # DV
            self._ref.shape[2] // 2,  # ML
        ]

    @property
    def ref(self) -> np.ndarray:
        return self._ref

    @property
    def ann(self) -> np.ndarray:
        return self._ann

    @property
    def plane(self) -> PlaneMode:
        return self._plane

    @property
    def n_ap_slices(self) -> int:
        return self._ref.shape[0]

    @property
    def n_dv_slices(self) -> int:
        return self._ref.shape[1]

    @property
    def n_ml_slices(self) -> int:
        return self._ref.shape[2]

    @property
    def shape_dv(self) -> int:
        return self._ref.shape[1]

    @property
    def shape_ml(self) -> int:
        return self._ref.shape[2]

    @property
    def ap_index(self) -> int:
        return self._slice_indices[0]

    @property
    def res_ap(self) -> float:
        return float(self.resolution[0])

    @property
    def res_dv(self) -> float:
        return float(self.resolution[1])

    @property
    def res_ml(self) -> float:
        return float(self.resolution[2])

    def set_ap_index(self, idx: int) -> None:
        self._slice_indices[0] = int(np.clip(idx, 0, self.n_ap_slices - 1))

    def set_slice_index(self, axis: int, idx: int) -> None:
        n = self._ref.shape[axis]
        self._slice_indices[axis] = int(np.clip(idx, 0, n - 1))

    def get_slice_index(self, axis: int) -> int:
        return self._slice_indices[axis]

    def n_slices_for_plane(self) -> int:
        """Number of slices for the current plane's slice axis."""
        cfg = PLANE_CONFIG[self._plane]
        return self._ref.shape[cfg["slice_axis"]]

    def current_slice_index(self) -> int:
        """Current slice index for the plane's slice axis."""
        cfg = PLANE_CONFIG[self._plane]
        return self._slice_indices[cfg["slice_axis"]]

    def set_current_slice_index(self, idx: int) -> None:
        cfg = PLANE_CONFIG[self._plane]
        self.set_slice_index(cfg["slice_axis"], idx)

    def get_displayed_resolution(self) -> tuple[float, float]:
        """(res_axis0, res_axis1) in μm for the displayed plane."""
        r = self.resolution
        return PLANE_CONFIG[self._plane]["res"](r)

    def reset_to_center(self) -> None:
        """Reset slice indices to center of each axis."""
        self._slice_indices = [
            self._ref.shape[0] // 2,
            self._ref.shape[1] // 2,
            self._ref.shape[2] // 2,
        ]

    def set_dims_range(self, viewer: Viewer) -> None:
        """Explicitly set dims range to atlas extent. Prevents incorrect range from 2D shapes layer."""
        for axis in range(3):
            n = self._ref.shape[axis]
            res = float(self.resolution[axis])
            # World coords: 0 to (n-1)*res, step res -> n steps (0..n-1)
            start = 0.0
            stop = max(start, (n - 1) * res)
            viewer.dims.set_range(axis, (start, stop, res))

    def apply_plane(self, viewer: Viewer, plane: PlaneMode) -> None:
        """Switch viewer to the given plane. Updates dims order, labels, current step."""
        self._plane = plane
        cfg = PLANE_CONFIG[plane]
        viewer.dims.order = cfg["order"]
        viewer.dims.axis_labels = ("AP", "DV", "ML")
        self.set_dims_range(viewer)
        step = tuple(self._slice_indices)
        viewer.dims.current_step = step
        viewer.dims.ndisplay = 2

    def sync_step_from_viewer(self, step: tuple[int, ...]) -> None:
        """Update our slice indices from viewer's current_step (AP, DV, ML order)."""
        for i, v in enumerate(step[:3]):
            self._slice_indices[i] = int(v)

    def add_to_viewer(self, viewer: Viewer) -> tuple[str, str]:
        """Add 3D reference and annotation as napari layers. Configure for coronal view.

        Returns (ref_layer_name, ann_layer_name).
        """
        ref_layer = viewer.add_image(
            self._ref,
            name=f"{self.atlas.atlas_name}_reference",
            scale=self.resolution,
            visible=True,
        )
        ann_layer = viewer.add_labels(
            self._ann,
            name=f"{self.atlas.atlas_name}_annotation",
            scale=self.resolution,
        )
        self.apply_plane(viewer, "coronal")
        return ref_layer.name, ann_layer.name

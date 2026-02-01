"""Atlas coronal view and layer management."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from brainglobe_atlasapi import BrainGlobeAtlas

if TYPE_CHECKING:
    from napari.viewer import Viewer


class AtlasCoronalView:
    """Manages Allen mouse brain atlas and coronal slice display.

    Atlas axes are (AP, DV, ML): axis 0 = AP, axis 1 = DV, axis 2 = ML.
    Coronal slice = constant AP, displays (DV, ML) plane.
    """

    def __init__(self, atlas_name: str) -> None:
        self.atlas = BrainGlobeAtlas(atlas_name=atlas_name)
        self._ref = np.asarray(self.atlas.reference)
        self._ann = np.asarray(self.atlas.annotation)
        self.resolution = self.atlas.resolution  # (res_ap, res_dv, res_ml) μm
        self.ap_index = self._ref.shape[0] // 2

    @property
    def ref(self) -> np.ndarray:
        return self._ref

    @property
    def ann(self) -> np.ndarray:
        return self._ann

    @property
    def n_ap_slices(self) -> int:
        return self._ref.shape[0]

    @property
    def shape_dv(self) -> int:
        return self._ref.shape[1]

    @property
    def shape_ml(self) -> int:
        return self._ref.shape[2]

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
        self.ap_index = int(np.clip(idx, 0, self.n_ap_slices - 1))

    def add_to_viewer(self, viewer: Viewer) -> tuple[str, str]:
        """Add 3D reference and annotation as napari layers. Configure for coronal view.

        Returns (ref_layer_name, ann_layer_name).
        """
        # Add 3D reference and annotation
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

        # Coronal view: display (DV, ML), slice along AP
        # napari dims: last dims are displayed. We want (DV, ML) displayed, AP as slider.
        # Default order for 3D is (0, 1, 2) -> display (1, 2) = (DV, ML), slice 0 = AP.
        viewer.dims.order = (0, 1, 2)
        viewer.dims.axis_labels = ("AP", "DV", "ML")
        viewer.dims.current_step = (self.ap_index, 0, 0)

        # Start in 2D mode for coronal viewing
        viewer.dims.ndisplay = 2

        return ref_layer.name, ann_layer.name

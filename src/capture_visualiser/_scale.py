"""μm ↔ voxel scaling utilities."""

from __future__ import annotations


def mm_to_voxels(
    width_mm: float,
    height_mm: float,
    res_dv_um: float,
    res_ml_um: float,
) -> tuple[float, float]:
    """Convert capture area dimensions from mm to voxels.

    width_mm maps to ML axis, height_mm maps to DV axis (coronal view).
    """
    width_um = width_mm * 1000.0
    height_um = height_mm * 1000.0
    width_vox = width_um / res_ml_um
    height_vox = height_um / res_dv_um
    return width_vox, height_vox

"""Capture area rectangle overlay for spatial transcriptomics platforms."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from capture_visualiser._scale import mm_to_voxels

if TYPE_CHECKING:
    from napari.viewer import Viewer

# Distinct colors for multiple capture areas
CAPTURE_COLORS = [
    "cyan",
    "magenta",
    "lime",
    "orange",
    "gold",
    "deepskyblue",
    "hotpink",
    "springgreen",
    "coral",
    "mediumpurple",
]


def create_capture_rectangle_vertices(
    center_dv: float,
    center_ml: float,
    width_vox: float,
    height_vox: float,
    rotation_rad: float = 0.0,
) -> np.ndarray:
    """Create rectangle vertices in (DV, ML) voxel coordinates.

    Rectangle is axis-aligned when rotation=0: width along ML, height along DV.
    Vertices order: bottom-left, bottom-right, top-right, top-left (closed polygon).
    """
    hx = width_vox / 2.0
    hy = height_vox / 2.0
    # Local coords: (ML, DV) for each vertex
    local = np.array(
        [
            [-hx, -hy],  # bottom-left
            [hx, -hy],   # bottom-right
            [hx, hy],    # top-right
            [-hx, hy],   # top-left
        ]
    )
    if rotation_rad != 0:
        c, s = math.cos(rotation_rad), math.sin(rotation_rad)
        R = np.array([[c, -s], [s, c]])
        local = local @ R.T
    # Convert to (DV, ML) for napari: swap columns
    # napari 2D: first axis = row (DV), second = col (ML)
    vertices = np.column_stack([local[:, 1] + center_dv, local[:, 0] + center_ml])
    return vertices


def add_capture_rectangle(
    viewer: Viewer,
    width_mm: float,
    height_mm: float,
    res_dv_um: float,
    res_ml_um: float,
    center_dv_vox: float | None = None,
    center_ml_vox: float | None = None,
    shape_dv: int | None = None,
    shape_ml: int | None = None,
    layer_name: str = "Capture areas",
    face_color: str = "cyan",
) -> str:
    """Add a Shapes layer with the capture rectangle, scaled to atlas.

    If center not given, places rectangle at center of (shape_dv, shape_ml).
    Returns the layer name.
    """
    width_vox, height_vox = mm_to_voxels(
        width_mm, height_mm, res_dv_um, res_ml_um
    )

    if center_dv_vox is None or center_ml_vox is None:
        if shape_dv is None or shape_ml is None:
            raise ValueError("Provide center or (shape_dv, shape_ml)")
        center_dv_vox = shape_dv / 2.0
        center_ml_vox = shape_ml / 2.0

    vertices = create_capture_rectangle_vertices(
        center_dv_vox, center_ml_vox, width_vox, height_vox, 0.0
    )
    # napari rectangle: 4 vertices as (N, 2) in data coords
    # scale must match atlas so rectangle overlays correctly
    # Scale matches atlas in-plane (DV, ML) so rectangle overlays correctly
    shapes_layer = viewer.add_shapes(
        [vertices],
        shape_type="polygon",
        scale=(res_dv_um, res_ml_um),
        name=layer_name,
        edge_color=face_color,
        face_color=face_color,
        opacity=0.3,
        edge_width=0,
    )
    return shapes_layer.name


def add_capture_shape_to_layer(
    shapes_layer,
    width_mm: float,
    height_mm: float,
    res_dv_um: float,
    res_ml_um: float,
    face_color: str,
    all_face_colors: list[str],
    center_dv_vox: float | None = None,
    center_ml_vox: float | None = None,
    shape_dv: int | None = None,
    shape_ml: int | None = None,
) -> None:
    """Add a new capture rectangle to an existing Shapes layer.

    all_face_colors must be the full list of colors for all shapes after adding
    (existing colors + face_color).
    """
    width_vox, height_vox = mm_to_voxels(
        width_mm, height_mm, res_dv_um, res_ml_um
    )
    if center_dv_vox is None or center_ml_vox is None:
        if shape_dv is None or shape_ml is None:
            raise ValueError("Provide center or (shape_dv, shape_ml)")
        center_dv_vox = shape_dv / 2.0
        center_ml_vox = shape_ml / 2.0

    vertices = create_capture_rectangle_vertices(
        center_dv_vox, center_ml_vox, width_vox, height_vox, 0.0
    )
    new_data = list(shapes_layer.data) + [vertices]
    shapes_layer.data = new_data
    shapes_layer.face_color = all_face_colors
    shapes_layer.edge_color = all_face_colors


def update_capture_shape_at_index(
    shapes_layer,
    shape_index: int,
    width_mm: float,
    height_mm: float,
    res_dv_um: float,
    res_ml_um: float,
) -> None:
    """Update the capture rectangle at shape_index, preserving its center and rotation."""
    if shapes_layer.nshapes == 0 or shape_index >= shapes_layer.nshapes:
        return
    vertices = np.asarray(shapes_layer.data[shape_index])
    # vertices are (DV, ML), 4 points
    center_dv = float(np.mean(vertices[:, 0]))
    center_ml = float(np.mean(vertices[:, 1]))
    # Rotation: angle of edge 0->1 from horizontal (ML axis)
    edge = vertices[1] - vertices[0]
    rotation_rad = math.atan2(edge[0], edge[1])
    width_vox, height_vox = mm_to_voxels(
        width_mm, height_mm, res_dv_um, res_ml_um
    )
    new_vertices = create_capture_rectangle_vertices(
        center_dv, center_ml, width_vox, height_vox, rotation_rad
    )
    new_data = list(shapes_layer.data)
    new_data[shape_index] = new_vertices
    shapes_layer.data = new_data


def rotate_capture_shape_at_index(
    shapes_layer,
    shape_index: int,
    angle_deg: float,
    res_dv_um: float,
    res_ml_um: float,
) -> None:
    """Rotate the capture rectangle at shape_index by angle_deg degrees around its center."""
    if shapes_layer.nshapes == 0 or shape_index >= shapes_layer.nshapes:
        return
    vertices = np.asarray(shapes_layer.data[shape_index])
    center_dv = float(np.mean(vertices[:, 0]))
    center_ml = float(np.mean(vertices[:, 1]))
    edge = vertices[1] - vertices[0]
    rotation_rad = math.atan2(edge[0], edge[1])
    # Current size from vertices: edge lengths in voxels
    w_vox = float(np.linalg.norm(vertices[1] - vertices[0]))
    h_vox = float(np.linalg.norm(vertices[3] - vertices[0]))
    new_rotation_rad = rotation_rad + math.radians(angle_deg)
    new_vertices = create_capture_rectangle_vertices(
        center_dv, center_ml, w_vox, h_vox, new_rotation_rad
    )
    new_data = list(shapes_layer.data)
    new_data[shape_index] = new_vertices
    shapes_layer.data = new_data


def get_shapes_physical_params(shapes_layer, res_dv_um: float, res_ml_um: float):
    """Extract shape parameters in physical units (μm) for atlas switching.

    Returns list of dicts: center_dv_um, center_ml_um, width_um, height_um, rotation_rad.
    """
    params_list = []
    for i in range(shapes_layer.nshapes):
        vertices = np.asarray(shapes_layer.data[i])
        center_dv = float(np.mean(vertices[:, 0]))
        center_ml = float(np.mean(vertices[:, 1]))
        center_dv_um = center_dv * res_dv_um
        center_ml_um = center_ml * res_ml_um
        edge_w = vertices[1] - vertices[0]
        edge_h = vertices[3] - vertices[0]
        width_um = float(np.sqrt((edge_w[0] * res_dv_um) ** 2 + (edge_w[1] * res_ml_um) ** 2))
        height_um = float(np.sqrt((edge_h[0] * res_dv_um) ** 2 + (edge_h[1] * res_ml_um) ** 2))
        rotation_rad = math.atan2(edge_w[0], edge_w[1])
        params_list.append({
            "center_dv_um": center_dv_um,
            "center_ml_um": center_ml_um,
            "width_um": width_um,
            "height_um": height_um,
            "rotation_rad": rotation_rad,
        })
    return params_list


def recreate_shapes_in_new_atlas(
    shapes_layer,
    params_list: list,
    res_dv_um: float,
    res_ml_um: float,
    colors: list[str],
) -> None:
    """Recreate shapes from physical params in new atlas coordinate space."""
    new_data = []
    for p in params_list:
        center_dv_vox = p["center_dv_um"] / res_dv_um
        center_ml_vox = p["center_ml_um"] / res_ml_um
        width_vox = p["width_um"] / res_ml_um
        height_vox = p["height_um"] / res_dv_um
        verts = create_capture_rectangle_vertices(
            center_dv_vox, center_ml_vox, width_vox, height_vox,
            p["rotation_rad"]
        )
        new_data.append(verts)
    shapes_layer.data = new_data
    shapes_layer.scale = (res_dv_um, res_ml_um)
    shapes_layer.face_color = colors
    shapes_layer.edge_color = colors
    shapes_layer.edge_width = 0

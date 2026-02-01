"""Capture area rectangle overlay for spatial transcriptomics platforms."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Literal

import numpy as np

from capture_visualiser._constants import XENIUM_SLIDE_HEIGHT_MM, XENIUM_SLIDE_WIDTH_MM
from capture_visualiser._scale import mm_to_voxels

if TYPE_CHECKING:
    from napari.viewer import Viewer

PlaneMode = Literal["coronal", "sagittal", "horizontal"]

# Default AP extent (μm) for capture area when viewed in sagittal/horizontal
DEFAULT_DEPTH_AP_UM = 1000.0  # 1 mm

# Shapes per capture area: 1 planar rectangle (visible in 2D slice view)
SHAPES_PER_CAPTURE = 1

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


def create_rectangle_vertices_2d(
    center_axis0: float,
    center_axis1: float,
    height_axis0: float,
    width_axis1: float,
    rotation_rad: float = 0.0,
) -> np.ndarray:
    """Create rectangle vertices in (axis0, axis1) coordinates.

    height along axis0, width along axis1. rotation_rad rotates in the plane.
    Returns (N, 2) array.
    """
    h0 = height_axis0 / 2.0
    h1 = width_axis1 / 2.0
    local = np.array(
        [
            [-h1, -h0],  # bottom-left
            [h1, -h0],   # bottom-right
            [h1, h0],    # top-right
            [-h1, h0],   # top-left
        ]
    )
    if rotation_rad != 0:
        c, s = math.cos(rotation_rad), math.sin(rotation_rad)
        R = np.array([[c, -s], [s, c]])
        local = local @ R.T
    vertices = np.column_stack([local[:, 1] + center_axis0, local[:, 0] + center_axis1])
    return vertices


def create_capture_rectangle_vertices(
    center_dv: float,
    center_ml: float,
    width_vox: float,
    height_vox: float,
    rotation_rad: float = 0.0,
) -> np.ndarray:
    """Create rectangle vertices in (DV, ML) voxel coordinates.

    Rectangle is axis-aligned when rotation=0: width along ML, height along DV.
    """
    return create_rectangle_vertices_2d(
        center_dv, center_ml, height_vox, width_vox, rotation_rad
    )


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


def add_capture_shape_to_layer_3d(
    shapes_layer,
    width_mm: float,
    height_mm: float,
    plane: PlaneMode,
    res_ap_um: float,
    res_dv_um: float,
    res_ml_um: float,
    face_color: str,
    all_face_colors: list[str],
    shape_ap: int,
    shape_dv: int,
    shape_ml: int,
    slice_index: int,
) -> None:
    """Add a new 3D capture rectangle (planar, in slice) to an existing Shapes layer."""
    if plane == "coronal":
        res_0, res_1 = res_dv_um, res_ml_um
        slice_res = res_ap_um
        center_0, center_1 = shape_dv / 2.0, shape_ml / 2.0
        w_mm, h_mm = width_mm, height_mm
    elif plane == "sagittal":
        res_0, res_1 = res_ap_um, res_dv_um
        slice_res = res_ml_um
        center_0, center_1 = shape_ap / 2.0, shape_dv / 2.0
        w_mm, h_mm = width_mm, height_mm
    else:
        res_0, res_1 = res_ap_um, res_ml_um
        slice_res = res_dv_um
        center_0, center_1 = shape_ap / 2.0, shape_ml / 2.0
        w_mm, h_mm = width_mm, height_mm
    width_vox = (w_mm * 1000) / res_1
    height_vox = (h_mm * 1000) / res_0
    verts = _create_planar_rectangle(
        plane, center_0, center_1, height_vox, width_vox,
        float(slice_index), 0.0,
    )
    new_data = list(shapes_layer.data) + [verts]
    shapes_layer.data = new_data
    shapes_layer.face_color = all_face_colors
    shapes_layer.edge_color = all_face_colors


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


def update_capture_shape_at_index_3d(
    shapes_layer,
    shape_indices: list[int],
    width_mm: float,
    height_mm: float,
    plane: PlaneMode,
    res_ap_um: float,
    res_dv_um: float,
    res_ml_um: float,
) -> None:
    """Update 3D shape(s) at indices (cuboid group), preserving center and rotation."""
    if not shape_indices or shapes_layer.nshapes == 0:
        return
    shape_index = shape_indices[0]
    if shape_index >= shapes_layer.nshapes:
        return
    vertices = np.asarray(shapes_layer.data[shape_index])
    if vertices.shape[1] < 3:
        for idx in shape_indices:
            update_capture_shape_at_index(
                shapes_layer, idx, width_mm, height_mm, res_dv_um, res_ml_um
            )
        return
    if plane == "coronal":
        disp_idx, res_0, res_1 = (1, 2), res_dv_um, res_ml_um
        slice_res = res_ap_um
        w_mm, h_mm = width_mm, height_mm
    elif plane == "sagittal":
        disp_idx, res_0, res_1 = (0, 1), res_ap_um, res_dv_um
        slice_res = res_ml_um
        w_mm, h_mm = width_mm, height_mm
    else:
        disp_idx, res_0, res_1 = (0, 2), res_ap_um, res_ml_um
        slice_res = res_dv_um
        w_mm, h_mm = width_mm, height_mm
    v0 = float(np.mean(vertices[:, disp_idx[0]]))
    v1 = float(np.mean(vertices[:, disp_idx[1]]))
    slice_idx = 0 if plane == "coronal" else 2 if plane == "sagittal" else 1
    slice_center = float(vertices[0, slice_idx])
    edge = vertices[1, disp_idx] - vertices[0, disp_idx]
    rotation_rad = math.atan2(edge[0], edge[1])
    width_vox = (w_mm * 1000) / res_1
    height_vox = (h_mm * 1000) / res_0
    verts = _create_planar_rectangle(
        plane, v0, v1, height_vox, width_vox,
        slice_center, rotation_rad,
    )
    new_data = list(shapes_layer.data)
    if shape_indices:
        new_data[shape_indices[0]] = verts
    shapes_layer.data = new_data


def rotate_capture_shape_at_index_3d(
    shapes_layer,
    shape_indices: list[int],
    angle_deg: float,
    plane: PlaneMode,
    res_ap_um: float,
    res_dv_um: float,
    res_ml_um: float,
) -> None:
    """Rotate 3D shape(s) at indices (cuboid group) by angle_deg in the displayed plane."""
    if not shape_indices or shapes_layer.nshapes == 0:
        return
    shape_index = shape_indices[0]
    if shape_index >= shapes_layer.nshapes:
        return
    vertices = np.asarray(shapes_layer.data[shape_index])
    if vertices.shape[1] < 3:
        for idx in shape_indices:
            rotate_capture_shape_at_index(
                shapes_layer, idx, angle_deg, res_dv_um, res_ml_um
            )
        return
    if plane == "coronal":
        disp_idx, res_0, res_1 = (1, 2), res_dv_um, res_ml_um
        slice_res = res_ap_um
    elif plane == "sagittal":
        disp_idx, res_0, res_1 = (0, 1), res_ap_um, res_dv_um
        slice_res = res_ml_um
    else:
        disp_idx, res_0, res_1 = (0, 2), res_ap_um, res_ml_um
        slice_res = res_dv_um
    v0 = float(np.mean(vertices[:, disp_idx[0]]))
    v1 = float(np.mean(vertices[:, disp_idx[1]]))
    slice_idx = 0 if plane == "coronal" else 2 if plane == "sagittal" else 1
    slice_center = float(vertices[0, slice_idx])
    edge = vertices[1, disp_idx] - vertices[0, disp_idx]
    rotation_rad = math.atan2(edge[0], edge[1])
    w_vox = float(np.linalg.norm(vertices[1, disp_idx] - vertices[0, disp_idx]))
    h_vox = float(np.linalg.norm(vertices[3, disp_idx] - vertices[0, disp_idx]))
    new_rotation_rad = rotation_rad + math.radians(angle_deg)
    verts = _create_planar_rectangle(
        plane, v0, v1, h_vox, w_vox,
        slice_center, new_rotation_rad,
    )
    new_data = list(shapes_layer.data)
    if shape_indices:
        new_data[shape_indices[0]] = verts
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


def get_shapes_physical_params(
    shapes_layer,
    plane: PlaneMode,
    res_ap_um: float,
    res_dv_um: float,
    res_ml_um: float,
    slice_position_um: float,
    prev_params: list[dict] | None = None,
) -> list[dict]:
    """Extract shape parameters in physical units (μm) for atlas/plane switching.

    slice_position_um: current slice position along slice axis (μm).
    prev_params: when switching from sagittal/horizontal, preserves width_ml/height_dv.
    Returns list of dicts with full 3D params.
    """
    params_list = []
    step = SHAPES_PER_CAPTURE
    for i in range(0, shapes_layer.nshapes, step):
        area_idx = i // step
        prev = prev_params[area_idx] if prev_params and area_idx < len(prev_params) else {}
        vertices = np.asarray(shapes_layer.data[i])
        is_3d = vertices.shape[1] >= 3
        if is_3d:
            if plane == "coronal":
                disp_idx = (1, 2)
                slice_idx, slice_res = 0, res_ap_um
            elif plane == "sagittal":
                disp_idx = (0, 1)
                slice_idx, slice_res = 2, res_ml_um
            else:
                disp_idx = (0, 2)
                slice_idx, slice_res = 1, res_dv_um
            v0 = float(np.mean(vertices[:, disp_idx[0]]))
            v1 = float(np.mean(vertices[:, disp_idx[1]]))
            slice_pos = float(vertices[0, slice_idx])
            slice_position_um = slice_pos * slice_res
            edge_w = vertices[1, disp_idx] - vertices[0, disp_idx]
            edge_h = vertices[3, disp_idx] - vertices[0, disp_idx]
        else:
            v0 = float(np.mean(vertices[:, 0]))
            v1 = float(np.mean(vertices[:, 1]))
            slice_position_um = slice_position_um  # use param
            edge_w = vertices[1] - vertices[0]
            edge_h = vertices[3] - vertices[0]

        if plane == "coronal":
            res_0, res_1 = res_dv_um, res_ml_um
            center_dv_um = v0 * res_0
            center_ml_um = v1 * res_1
            center_ap_um = slice_position_um
            width_ml_um = float(np.sqrt((edge_w[0] * res_0) ** 2 + (edge_w[1] * res_1) ** 2))
            height_dv_um = float(np.sqrt((edge_h[0] * res_0) ** 2 + (edge_h[1] * res_1) ** 2))
            depth_ap_um = prev.get("depth_ap_um", DEFAULT_DEPTH_AP_UM)
            rotation_rad = math.atan2(edge_w[0], edge_w[1])
        elif plane == "sagittal":
            res_0, res_1 = res_ap_um, res_dv_um
            center_ap_um = v0 * res_0
            center_dv_um = v1 * res_1
            center_ml_um = slice_position_um
            depth_ap_um = float(np.sqrt((edge_w[0] * res_0) ** 2 + (edge_w[1] * res_1) ** 2))
            height_dv_um = float(np.sqrt((edge_h[0] * res_0) ** 2 + (edge_h[1] * res_1) ** 2))
            width_ml_um = prev.get("width_ml_um", height_dv_um)
            rotation_rad = prev.get("rotation_rad", 0.0)
        else:  # horizontal
            res_0, res_1 = res_ap_um, res_ml_um
            center_ap_um = v0 * res_0
            center_ml_um = v1 * res_1
            center_dv_um = slice_position_um
            depth_ap_um = float(np.sqrt((edge_w[0] * res_0) ** 2 + (edge_w[1] * res_1) ** 2))
            width_ml_um = float(np.sqrt((edge_h[0] * res_0) ** 2 + (edge_h[1] * res_1) ** 2))
            height_dv_um = prev.get("height_dv_um", depth_ap_um)
            rotation_rad = prev.get("rotation_rad", 0.0)

        p = {
            "center_ap_um": center_ap_um,
            "center_dv_um": center_dv_um,
            "center_ml_um": center_ml_um,
            "width_ml_um": width_ml_um,
            "height_dv_um": height_dv_um,
            "depth_ap_um": depth_ap_um,
            "rotation_rad": rotation_rad,
        }
        params_list.append(p)
    return params_list


def recreate_shapes_for_plane(
    shapes_layer,
    params_list: list[dict],
    plane: PlaneMode,
    res_ap_um: float,
    res_dv_um: float,
    res_ml_um: float,
    colors: list[str],
    use_3d: bool = True,
) -> None:
    """Recreate shapes in the given plane from 3D physical params. use_3d=True for 3D vertices."""
    if plane == "coronal":
        res_0, res_1 = res_dv_um, res_ml_um
    elif plane == "sagittal":
        res_0, res_1 = res_ap_um, res_dv_um
    else:
        res_0, res_1 = res_ap_um, res_ml_um

    new_data = []
    new_colors = []
    for p_idx, p in enumerate(params_list):
        cap = p["center_ap_um"]
        cdv = p["center_dv_um"]
        cml = p["center_ml_um"]
        w_ml = p["width_ml_um"]
        h_dv = p["height_dv_um"]
        d_ap = p.get("depth_ap_um", DEFAULT_DEPTH_AP_UM)
        rot = p["rotation_rad"]

        if plane == "coronal":
            center_0, center_1 = cdv, cml
            size_0, size_1 = h_dv, w_ml
            slice_vox = cap / res_ap_um
            slice_res = res_ap_um
            rot_rad = rot
        elif plane == "sagittal":
            center_0, center_1 = cap, cdv
            size_0, size_1 = d_ap, h_dv
            slice_vox = cml / res_ml_um
            slice_res = res_ml_um
            rot_rad = 0.0
        else:
            center_0, center_1 = cap, cml
            size_0, size_1 = d_ap, w_ml
            slice_vox = cdv / res_dv_um
            slice_res = res_dv_um
            rot_rad = 0.0

        c0_vox = center_0 / res_0
        c1_vox = center_1 / res_1
        s0_vox = size_0 / res_0
        s1_vox = size_1 / res_1
        if use_3d:
            verts = _create_planar_rectangle(
                plane, c0_vox, c1_vox, s0_vox, s1_vox,
                slice_vox, rot_rad,
            )
            new_data.append(verts)
            col = colors[p_idx] if p_idx < len(colors) else CAPTURE_COLORS[0]
            new_colors.append(col)
        else:
            verts_2d = create_rectangle_vertices_2d(c0_vox, c1_vox, s0_vox, s1_vox, rot_rad)
            new_data.append(verts_2d)
            new_colors.append(colors[p_idx] if p_idx < len(colors) else CAPTURE_COLORS[0])

    shapes_layer.data = new_data
    if use_3d:
        shapes_layer.scale = (res_ap_um, res_dv_um, res_ml_um)
        if new_colors:
            shapes_layer.face_color = new_colors
            shapes_layer.edge_color = new_colors
    else:
        shapes_layer.scale = (res_0, res_1)
        shapes_layer.face_color = colors
        shapes_layer.edge_color = colors
    shapes_layer.edge_width = 0


def update_capture_slice_positions(
    shapes_layer,
    plane: PlaneMode,
    current_slice_index: float,
) -> None:
    """Translate capture areas (matching plane only) so their slice-dimension = current_slice_index."""
    if shapes_layer.nshapes == 0:
        return
    slice_axis = 0 if plane == "coronal" else 2 if plane == "sagittal" else 1
    new_data = list(shapes_layer.data)
    for idx in range(shapes_layer.nshapes):
        v = np.asarray(new_data[idx])
        if v.shape[1] < 3:
            continue
        shape_plane = infer_shape_face_plane(v)
        if shape_plane != plane:
            continue
        current_slice_val = float(np.mean(v[:, slice_axis]))
        delta = current_slice_index - current_slice_val
        if abs(delta) < 0.01:
            continue
        v_new = v.copy()
        v_new[:, slice_axis] += delta
        new_data[idx] = v_new
    shapes_layer.data = new_data


def infer_shape_face_plane(vertices: np.ndarray) -> PlaneMode:
    """Infer which plane a 3D shape faces from its vertices.

    The dimension with smallest variance is the slice (constant) dimension.
    Axis 0=AP -> coronal, 1=DV -> horizontal, 2=ML -> sagittal.
    """
    verts = np.asarray(vertices)
    if verts.shape[0] < 2 or verts.shape[1] < 3:
        return "coronal"
    variances = np.var(verts[:, :3], axis=0)
    slice_axis = int(np.argmin(variances))
    return {0: "coronal", 1: "horizontal", 2: "sagittal"}[slice_axis]


def _verts_2d_to_3d(
    verts_2d: np.ndarray,
    plane: PlaneMode,
    slice_index: float,
) -> np.ndarray:
    """Convert 2D (axis0, axis1) vertices to 3D (AP, DV, ML)."""
    verts_3d = np.zeros((len(verts_2d), 3))
    if plane == "coronal":
        verts_3d[:, 0] = slice_index
        verts_3d[:, 1] = verts_2d[:, 0]  # DV
        verts_3d[:, 2] = verts_2d[:, 1]  # ML
    elif plane == "sagittal":
        verts_3d[:, 0] = verts_2d[:, 0]  # AP
        verts_3d[:, 1] = verts_2d[:, 1]  # DV
        verts_3d[:, 2] = slice_index
    else:  # horizontal
        verts_3d[:, 0] = verts_2d[:, 0]  # AP
        verts_3d[:, 1] = slice_index
        verts_3d[:, 2] = verts_2d[:, 1]  # ML
    return verts_3d


def _create_planar_rectangle(
    plane: PlaneMode,
    center_0: float,
    center_1: float,
    size_0_vox: float,
    size_1_vox: float,
    slice_center_vox: float,
    rotation_rad: float = 0.0,
) -> np.ndarray:
    """Create a single 4-vertex polygon in the slice plane (visible in 2D view)."""
    verts_2d = create_rectangle_vertices_2d(
        center_0, center_1, size_0_vox, size_1_vox, rotation_rad
    )
    return _verts_2d_to_3d(verts_2d, plane, slice_center_vox)


def add_default_capture_for_plane(
    viewer,
    plane: PlaneMode,
    res_ap_um: float,
    res_dv_um: float,
    res_ml_um: float,
    shape_ap: int,
    shape_dv: int,
    shape_ml: int,
    slice_index: int,
    layer_name: str = "Capture areas",
) -> str:
    """Add a new 3D Shapes layer with default Xenium capture (cuboid)."""
    if plane == "coronal":
        res_0, res_1 = res_dv_um, res_ml_um
        slice_res = res_ap_um
        center_0, center_1 = shape_dv / 2.0, shape_ml / 2.0
        size_0 = XENIUM_SLIDE_HEIGHT_MM * 1000
        size_1 = XENIUM_SLIDE_WIDTH_MM * 1000
    elif plane == "sagittal":
        res_0, res_1 = res_ap_um, res_dv_um
        slice_res = res_ml_um
        center_0, center_1 = shape_ap / 2.0, shape_dv / 2.0
        size_0 = DEFAULT_DEPTH_AP_UM
        size_1 = XENIUM_SLIDE_HEIGHT_MM * 1000
    else:
        res_0, res_1 = res_ap_um, res_ml_um
        slice_res = res_dv_um
        center_0, center_1 = shape_ap / 2.0, shape_ml / 2.0
        size_0 = DEFAULT_DEPTH_AP_UM
        size_1 = XENIUM_SLIDE_WIDTH_MM * 1000
    s0_vox = size_0 / res_0
    s1_vox = size_1 / res_1
    verts = _create_planar_rectangle(
        plane, center_0, center_1, s0_vox, s1_vox,
        float(slice_index), 0.0,
    )
    layer = viewer.add_shapes(
        [verts],
        shape_type="polygon",
        scale=(res_ap_um, res_dv_um, res_ml_um),
        ndim=3,
        name=layer_name,
        edge_color=CAPTURE_COLORS[0],
        face_color=CAPTURE_COLORS[0],
        opacity=0.3,
        edge_width=0,
    )
    return layer.name


def recreate_shapes_in_new_atlas(
    shapes_layer,
    params_list: list,
    plane: str,
    res_ap_um: float,
    res_dv_um: float,
    res_ml_um: float,
    colors: list[str],
) -> None:
    """Recreate 3D shapes from physical params in new atlas."""
    recreate_shapes_for_plane(
        shapes_layer,
        params_list,
        plane,
        res_ap_um,
        res_dv_um,
        res_ml_um,
        colors,
        use_3d=True,
    )

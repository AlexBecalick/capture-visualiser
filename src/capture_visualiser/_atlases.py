"""Atlas listing from brainglobe-atlasapi."""

from __future__ import annotations

# Recommended atlases (display "(recommended)" in UI)
RECOMMENDED_ATLASES = frozenset({"allen_mouse_25um", "allen_human_500um"})


def get_available_atlases() -> list[tuple[str, str]]:
    """Return list of (atlas_id, display_name) for all available atlases.

    display_name includes ' (recommended)' for recommended atlases.
    Sorted with recommended first, then alphabetically by id.
    """
    try:
        from brainglobe_atlasapi.list_atlases import get_all_atlases_lastversions

        atlases = get_all_atlases_lastversions()
    except Exception:
        return _fallback_atlas_list()

    if not atlases:
        return _fallback_atlas_list()

    items: list[tuple[str, str]] = []
    for atlas_id in sorted(atlases.keys()):
        suffix = " (recommended)" if atlas_id in RECOMMENDED_ATLASES else ""
        display = f"{atlas_id}{suffix}"
        items.append((atlas_id, display))

    # Sort: recommended first, then alphabetically
    def sort_key(item: tuple[str, str]) -> tuple[int, str]:
        atlas_id = item[0]
        is_rec = 0 if atlas_id in RECOMMENDED_ATLASES else 1
        return (is_rec, atlas_id)

    items.sort(key=sort_key)
    return items


def _fallback_atlas_list() -> list[tuple[str, str]]:
    """Fallback when brainglobe API is unavailable."""
    return [
        ("allen_mouse_25um", "allen_mouse_25um (recommended)"),
        ("allen_human_500um", "allen_human_500um (recommended)"),
        ("allen_mouse_10um", "allen_mouse_10um"),
        ("allen_mouse_50um", "allen_mouse_50um"),
        ("allen_mouse_100um", "allen_mouse_100um"),
    ]

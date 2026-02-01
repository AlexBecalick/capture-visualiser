"""Atlas loading with optional download progress."""

from __future__ import annotations

from qtpy.QtCore import QObject, QThread, Signal


def is_atlas_downloaded(atlas_name: str) -> bool:
    """Check if atlas is already downloaded."""
    from brainglobe_atlasapi.list_atlases import get_downloaded_atlases
    return atlas_name in get_downloaded_atlases()


class AtlasInstallWorker(QObject):
    """Worker that runs install_atlas in a thread, emitting progress."""

    progress = Signal(int, int)  # completed, total
    finished = Signal()
    error = Signal(str)

    def __init__(self, atlas_name: str) -> None:
        super().__init__()
        self._atlas_name = atlas_name

    def run(self) -> None:
        try:
            from brainglobe_atlasapi.update_atlases import install_atlas

            def on_progress(completed: int, total: int) -> None:
                self.progress.emit(completed, total)

            install_atlas(self._atlas_name, fn_update=on_progress)
            self.finished.emit()
        except Exception as e:
            self.error.emit(str(e))

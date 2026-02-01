"""AP plane slider widget."""

from __future__ import annotations

from qtpy.QtWidgets import QHBoxLayout, QLabel, QSlider, QWidget
from qtpy.QtCore import Qt, Signal


class APSliderWidget(QWidget):
    """Slider to step through anterior-posterior (AP) coronal planes."""

    valueChanged = Signal(int)

    def __init__(
        self,
        min_val: int = 0,
        max_val: int = 100,
        value: int | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._min = min_val
        self._max = max_val
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self._label = QLabel("AP plane:")
        layout.addWidget(self._label)

        self._slider = QSlider(Qt.Orientation.Horizontal)
        self._slider.setMinimum(min_val)
        self._slider.setMaximum(max_val)
        if value is not None:
            self._slider.setValue(value)
        else:
            self._slider.setValue((min_val + max_val) // 2)
        self._slider.valueChanged.connect(self._on_slider_changed)
        layout.addWidget(self._slider, stretch=1)

        self._value_label = QLabel(str(self._slider.value()))
        self._value_label.setMinimumWidth(40)
        layout.addWidget(self._value_label)

    def _on_slider_changed(self, val: int) -> None:
        self._value_label.setText(str(val))
        self.valueChanged.emit(val)

    def set_range(self, min_val: int, max_val: int) -> None:
        self._min = min_val
        self._max = max_val
        self._slider.setMinimum(min_val)
        self._slider.setMaximum(max_val)

    def set_value(self, value: int) -> None:
        self._slider.blockSignals(True)
        self._slider.setValue(int(max(self._min, min(self._max, value))))
        self._value_label.setText(str(self._slider.value()))
        self._slider.blockSignals(False)

    def value(self) -> int:
        return self._slider.value()

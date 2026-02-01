"""Capture Visualiser: Overlay spatial transcriptomics capture areas onto brain atlas sections."""

import warnings

# Suppress known napari shapes layer warning (NumPy 'where' with 'out' in triangulation)
warnings.filterwarnings(
    "ignore",
    message=".*'where' used without 'out'.*",
    category=UserWarning,
)

# Capture Visualiser

A napari plugin for overlaying spatial transcriptomics capture areas onto brain sections. Compare platforms such as Xenium, Visium HD, MERSCOPE, and GeoMx by browsing Allen mouse brain atlas coronal sections and positioning scaled rectangles to see which brain regions fit within each capture area.

## Features

- **BrainGlobe atlases**: Load any atlas from brainglobe-atlasapi (mouse, human, rat, zebrafish, etc.)
- **View planes**: Switch between coronal (default), sagittal, and horizontal views
- **Slice slider**: Step through sections (AP, ML, or DV depending on view plane)
- **Capture area overlays**: Rectangles scaled to platform dimensions (Xenium, Visium HD, MERSCOPE, GeoMx, BARseq)
- **Interactive placement**: Select, translate and rotate capture areas using napari's shape tools
- **Multiple capture areas**: Add several rectangles with different colours; select one to resize with presets

## Installation

### Option 1: Conda environment (recommended)

```bash
# From the project directory
cd capture_visualiser
conda env create -f environment.yml
conda activate capture-visualiser
```

This creates a conda environment and runs `pip install -e .` to install the package and its dependencies from `pyproject.toml`.

If you already have the environment and see "No Qt bindings could be found", install PyQt:

```bash
conda activate capture-visualiser
conda install -c conda-forge pyqt
```

### Option 2: pip (existing environment)

```bash
pip install -e .
```

### Requirements

- Python ≥ 3.11
- napari ≥ 0.6.1
- brainglobe-atlasapi ≥ 2.2.0

The Allen atlas is downloaded automatically on first use. The conda environment includes PyQt for napari's GUI.

## Usage

### Launching the plugin

1. Start napari:
   ```bash
   napari
   ```

2. Open the Capture Visualiser widget:
   - **Plugins** → **Capture Visualiser** → **Capture Visualiser**
   - Or use the napari Plugins menu and search for "Capture Visualiser"

### Workflow

1. **Load atlas**  
   Select an atlas from the dropdown (default: allen_mouse_25um, recommended for mouse) and click **Load atlas**. Atlases marked "(recommended)" are suggested starting points. A download progress bar appears when installing a new atlas. You can switch to another atlas at any time—capture areas are preserved and repositioned in the new atlas space.

2. **Choose view plane**  
   Click **Coronal** (default), **Sagittal**, or **Horizontal** to switch views. The slice slider updates to browse the corresponding axis (AP, ML, or DV). Capture areas are planar rectangles visible in 2D view. Only capture areas created in the current view plane follow the slice as you scroll; others stay fixed. Capture areas are kept when switching planes and retain the orientation they had when created.

3. **Browse sections**  
   Use the slice slider to move through sections and find the region you care about.

4. **Use presets**  
   Click a preset button (Xenium, Visium HD, MERSCOPE, GeoMx, BARseq) to resize the *selected* capture area. Hover for a short description. Select a shape first with the **Select shapes** tool.

5. **Custom capture area**  
   Enter width and height in mm, then click **Apply** to resize the selected capture area.

6. **Multiple capture areas**  
   Click **Create a new capture area** to add another rectangle. Each new area gets a different colour. Select any shape to resize it with presets or Apply.

7. **Position the capture rectangles**  
   Use napari's **Select shapes** tool to move and rotate rectangles over the atlas. Select one to make it the target for preset/Apply resize.

8. **Rotate**  
   Click **Rotate current capture area 90°** to rotate the selected rectangle.

### Tips

- **Create** and **Rotate** work in any view plane (coronal, sagittal, horizontal).
- Keep the viewer in **2D** mode for section browsing.
- Use napari's built-in dimension slider in addition to the plugin slider if you prefer.
- You can toggle layer visibility (reference vs annotation) in the napari layer list.

## Configuration

| Setting           | Description                         | Default |
|-------------------|-------------------------------------|---------|
| Width             | Capture area width in mm            | 10.45 (Xenium) |
| Height            | Capture area height in mm           | 22.45 (Xenium) |
| Atlas             | Dropdown selection from brainglobe-atlasapi | allen_mouse_25um (default) |
| Presets           | Config file                         | `presets.json` in package; override at `~/.config/capture_visualiser/presets.json` |

## License

BSD-3-Clause

## Acknowledgements

- [BrainGlobe](https://brainglobe.info/) for brainglobe-atlasapi and the Allen atlas integration
- [brainrender-napari](https://github.com/brainglobe/brainrender-napari) for atlas visualisation patterns in napari
- [napari](https://napari.org/) for the viewer and plugin system

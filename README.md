# LUMIN

**Library-Utilized Mixture Quantification in NMR**

LUMIN is a library-based workflow for efficient 1D 1H NMR qualitative and quantitative analysis of complex mixtures.

This project is developed by the research group of **Prof. Haibin Qu**, College of Pharmaceutical Sciences, Zhejiang University.

## How to Use

### Option A: Standard quantification workflow

Run:

```bash
python task05_Standard.py
```

This mode performs routine compound quantification using the current template library.

### Option B: Full-search quantification workflow

Run:

```bash
python task05_FullSearch.py
```

This mode performs a broader search strategy and then quantifies matched compounds.

### Option C: Interactive threshold exploration

Run:

```bash
python Window_thres_estimate.py
```

This opens a GUI tool to help inspect and optimize cluster/window threshold behavior.

## Typical Workflow

1. Put sample spectrum files (`.csv`) into the `sample` folder.
2. Make sure compound templates (`.csv`) are available in the `templates` folder.
3. Run `task05_Standard.py` (or `task05_FullSearch.py`).
4. Check generated Excel quantification tables and Lumin-style visualization outputs.

## About Output

- **Excel file** (`0624quantified-*.xlsx`):
  - Compound-level quantification results
  - Peak positions, cluster info, absolute/normalized/relative values
- **HTML file** (`fitting_html/visualization-*.html`):
  - Interactive Lumin-style spectrum visualization
- **JSON file** (`fitting_html/visualization-*.json`):
  - Serialized visualization data for downstream use

## Folder Notes

- `templates/`
  - Template spectra library used for searching and fitting.
  - Keep folder name unchanged unless you also update script paths.
- `sample/`
  - Input sample spectra (`shift, intensity` format).
- `fitting_visualization/`
  - Utility scripts for loading and visualizing saved fitting results.

## Template Extension

This project is provided as an example workflow.

If you need to analyze specific compounds that are not included in the current library, please prepare corresponding template files and place them in the `templates` folder before running analysis.

## Environment

Install dependencies:

```bash
pip install -r requirements.txt
```

## Contact

For collaboration, licensing, or technical questions, please contact Prof. Haibin Qu's research group at Zhejiang University.

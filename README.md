# Adaptable heatmap builder

`heatmap.py` builds heatmaps for lnDEF-like experiment tables without editing Python code.
It generates heatmaps from H2DCFDA-style analysis and supports different doses, concentrations, sample names, cell lines, pH values and material families through a JSON config or the Colab notebook.

## Quick start

Run the current GdF3/GdEuF3/BiCe example:

```bash
python heatmap.py
```

Run with an editable config:

```bash
python heatmap.py --config heatmap_config_example.json
```

Create a fresh config template:

```bash
python heatmap.py --write-example-config my_heatmap_config.json
```

Preview the normalized dataset before plotting:

```bash
python heatmap.py --config my_heatmap_config.json --preview
```

## Colab workflow

Open `heatmap_colab.ipynb` in Google Colab, then upload:

- `heatmap.py`
- your CSV/XLSX experiment files
- optional `heatmap_config_example.json` or your edited JSON config

The notebook has form fields for common cases, so users can change file names, pH values, grouping rules, labels and output names without touching the plotting code.

## Supported input formats

### 1. Long/tidy table

One row per measurement:

```csv
Material,Family,pH,Conc_mg_mL,Dose_Gy,lnDEF
GdF3,GdF3 / GdEuF3,6.5,0.01,3,0.18
GdF3,GdF3 / GdEuF3,6.5,0.05,3,0.53
```

Common aliases such as `sample`, `concentration`, `dose`, `value` and `mean` are detected automatically.

### 2. Wide dose table

One row per material/concentration/pH, with doses as columns:

```csv
Material,pH,Conc_mg_mL,3_Gy,6_Gy,9_Gy
GdF3,6.5,0.01,0.18,-0.22,-0.55
```

### 3. Excel triplet table

The current Excel layout is supported:

- first column: dose
- then repeated groups of 3 columns: mean, SD, N
- concentration is written in the sample header, for example `BiCe 5% (0.005)`

For this format the config must provide `ph`, because each file usually represents one pH.

## Config fields that matter most

- `sources`: list of input files and formats.
- `family_rules`: maps sample names to plot groups, for example `"contains": "BiCe"`.
- `include_families` / `exclude_families`: keeps one source from mixing unwanted groups.
- `plot.output_dir`: where PNG files are saved.
- `plot.group_titles`: nice titles for plot groups.
- `plot.material_order`: custom order for rows in each plot group.
- `plot.panel_column`: default is `pH`; can be changed to `cell_line`, `condition` or another metadata column if the table contains it.

from __future__ import annotations

import argparse
import json
import math
import os
import re
import tempfile
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize, TwoSlopeNorm


SUB3 = "3"

DEFAULT_CONFIG: dict[str, Any] = {
    "sources": [
        {
            "path": "Protons BiCe GdEuF3 6_5.xlsx",
            "format": "triplet_excel",
            "ph": 6.5,
            "include_families": ["BiCe"],
        },
        {
            "path": "Protons BiCe GdEuF3 7_4.xlsx",
            "format": "triplet_excel",
            "ph": 7.4,
            "include_families": ["BiCe"],
        },
        {
            "path": "gdf3_gdeuf3_manual.csv",
            "format": "wide_dose",
            "family": "GdF3 / GdEuF3",
        },
    ],
    "family_rules": [
        {"contains": "BiCe", "family": "BiCe"},
        {"contains": "GdEu", "family": "GdF3 / GdEuF3"},
        {"contains": "GdF3", "family": "GdF3 / GdEuF3"},
    ],
    "default_family": "All samples",
    "plot": {
        "output_dir": ".",
        "value_label": r"$\ln(\mathrm{DEF}_{\mathrm{ROS}})$",
        "colorbar_note": "blue: ROS decrease vs control; red: ROS increase vs control",
        "cmap": "bwr",
        "center": 0,
        "global_color_scale": True,
        "annotate": True,
        "dose_label": "Dose, Gy",
        "concentration_label": "Conc, mg/mL",
        "panel_label": "pH",
        "dpi": 300,
        "group_titles": {
            "GdF3 / GdEuF3": "GdF3 / GdEuF3",
        },
        "material_order": {
            "GdF3 / GdEuF3": ["GdEu 20% F3", "GdEu 5% F3", "GdF3"],
        },
        "output_names": {
            "BiCe": "lnDEF_heatmaps_BiCe_pH65_vs_74.png",
            "GdF3 / GdEuF3": [
                "lnDEF_heatmaps_GdF3_GdEuF3_pH65_vs_74.png",
                "lnDEF_heatmaps_Xray_GdF3_GdEuF3_pH65_vs_74_title_only_clean.png",
            ],
        },
    },
}

COLUMN_ALIASES = {
    "material": ["material", "sample", "sample_name", "np", "nanoparticle", "particle"],
    "family": ["family", "group", "series", "line", "cell_line", "cellline"],
    "pH": ["pH", "ph"],
    "dose": ["dose", "dose_gy", "dose,gy", "dose gy", "radiation_dose"],
    "conc": [
        "conc",
        "concentration",
        "conc_mg_ml",
        "conc_mg_mL",
        "concentration_mg_ml",
        "mg_ml",
    ],
    "lnDEF": ["lnDEF", "ln_def", "ln(def)", "value", "mean", "effect", "response"],
}

MISSING_VALUES = {"#NUM!", "#DIV/0!", "#VALUE!", "", "NA", "N/A", "nan", "NaN"}


def to_float(value: Any) -> float | None:
    """Convert Excel/CSV values to float while tolerating blanks and comma decimals."""
    try:
        if value is None:
            return None
        if isinstance(value, str):
            value = value.strip().replace(",", ".")
            if value in MISSING_VALUES:
                return None
        if pd.isna(value):
            return None
        result = float(value)
        if math.isnan(result):
            return None
        return result
    except Exception:
        return None


def compact_key(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def clean_label(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value).strip())


def clean_excel_duplicate_suffix(value: Any) -> str:
    return re.sub(r"\.\d+$", "", clean_label(value))


def natural_key(value: Any) -> list[Any]:
    parts = re.split(r"(\d+(?:[.,]\d+)?)", str(value))
    key: list[Any] = []
    for part in parts:
        number = to_float(part)
        key.append(number if number is not None else part.lower())
    return key


def format_material_label(name: str) -> str:
    return re.sub(r"F\s*3\b", f"F{SUB3}", str(name))


def slugify(value: Any) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    return slug.strip("_") or "heatmap"


def read_table(path: str | Path, sheet_name: str | int | None = None) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".xls", ".xlsx", ".xlsm"}:
        return pd.read_excel(path, sheet_name=0 if sheet_name is None else sheet_name)
    if suffix in {".tsv", ".tab"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path, sep=None, engine="python")


def find_column(
    df: pd.DataFrame,
    canonical_name: str,
    requested_name: str | None = None,
    required: bool = True,
) -> str | None:
    if requested_name and requested_name in df.columns:
        return requested_name

    aliases = [canonical_name, *COLUMN_ALIASES.get(canonical_name, [])]
    available = {compact_key(col): col for col in df.columns}
    for alias in aliases:
        hit = available.get(compact_key(alias))
        if hit is not None:
            return hit

    if required:
        raise ValueError(
            f"Column for '{canonical_name}' was not found. Available columns: {list(df.columns)}"
        )
    return None


def extract_dose_from_column(name: Any) -> float | None:
    text = clean_excel_duplicate_suffix(name)
    match = re.search(r"(^|\s|_)(\d+(?:[.,]\d+)?)(\s|_)*(gy|гр)?$", text, re.IGNORECASE)
    if not match:
        return None
    return to_float(match.group(2))


def infer_family(material: str, source: dict[str, Any], config: dict[str, Any]) -> str:
    if source.get("family"):
        return clean_label(source["family"])

    for rule in config.get("family_rules", []):
        family = rule.get("family")
        if not family:
            continue
        if "regex" in rule and re.search(str(rule["regex"]), material, flags=re.IGNORECASE):
            return clean_label(family)
        if "startswith" in rule and material.lower().startswith(str(rule["startswith"]).lower()):
            return clean_label(family)
        if "contains" in rule and str(rule["contains"]).lower() in material.lower():
            return clean_label(family)

    return clean_label(config.get("default_family", "All samples"))


def add_optional_metadata(record: dict[str, Any], row: pd.Series | None, source: dict[str, Any]) -> None:
    for key in ("experiment", "cell_line", "condition"):
        if source.get(key) is not None:
            record[key] = clean_label(source[key])
        elif row is not None:
            column = source.get(f"{key}_column")
            if column and column in row:
                record[key] = clean_label(row[column])


def load_long_table(source: dict[str, Any], config: dict[str, Any]) -> pd.DataFrame:
    df = read_table(source["path"], source.get("sheet_name"))
    columns = source.get("columns", {})

    material_col = find_column(df, "material", columns.get("material"))
    ph_col = find_column(df, "pH", columns.get("pH"), required=source.get("ph") is None)
    dose_col = find_column(df, "dose", columns.get("dose"))
    conc_col = find_column(df, "conc", columns.get("conc"))
    value_col = find_column(df, "lnDEF", columns.get("value") or columns.get("lnDEF"))
    family_col = find_column(df, "family", columns.get("family"), required=False)

    records: list[dict[str, Any]] = []
    for _, row in df.iterrows():
        material = clean_label(row[material_col])
        value = to_float(row[value_col])
        dose = to_float(row[dose_col])
        conc = to_float(row[conc_col])
        ph = to_float(source.get("ph")) if source.get("ph") is not None else to_float(row[ph_col])
        if value is None or dose is None or conc is None or ph is None:
            continue

        family = clean_label(row[family_col]) if family_col else infer_family(material, source, config)
        record = {
            "material": material,
            "family": family,
            "pH": ph,
            "dose": dose,
            "conc": conc,
            "lnDEF": value,
        }
        add_optional_metadata(record, row, source)
        records.append(record)

    return pd.DataFrame(records)


def load_wide_dose_table(source: dict[str, Any], config: dict[str, Any]) -> pd.DataFrame:
    df = read_table(source["path"], source.get("sheet_name"))
    columns = source.get("columns", {})

    material_col = find_column(df, "material", columns.get("material"))
    ph_col = find_column(df, "pH", columns.get("pH"), required=source.get("ph") is None)
    conc_col = find_column(df, "conc", columns.get("conc"))
    family_col = find_column(df, "family", columns.get("family"), required=False)

    if source.get("dose_columns"):
        dose_columns = {
            column_name: to_float(dose_value)
            for column_name, dose_value in source["dose_columns"].items()
        }
    else:
        dose_columns = {
            column_name: extract_dose_from_column(column_name)
            for column_name in df.columns
            if extract_dose_from_column(column_name) is not None
        }

    if not dose_columns:
        raise ValueError(
            f"No dose columns were detected in {source['path']}. "
            "Use columns like 3_Gy, 6_Gy or provide source['dose_columns']."
        )

    records: list[dict[str, Any]] = []
    for _, row in df.iterrows():
        material = clean_label(row[material_col])
        ph = to_float(source.get("ph")) if source.get("ph") is not None else to_float(row[ph_col])
        conc = to_float(row[conc_col])
        if ph is None or conc is None:
            continue
        family = clean_label(row[family_col]) if family_col else infer_family(material, source, config)

        for dose_col, dose in dose_columns.items():
            value = to_float(row[dose_col])
            if value is None or dose is None:
                continue
            record = {
                "material": material,
                "family": family,
                "pH": ph,
                "dose": dose,
                "conc": conc,
                "lnDEF": value,
            }
            add_optional_metadata(record, row, source)
            records.append(record)

    return pd.DataFrame(records)


def parse_concentration_from_label(label: str) -> float | None:
    match = re.search(r"\(([^)]+)\)", label)
    return to_float(match.group(1)) if match else None


def load_triplet_excel(source: dict[str, Any], config: dict[str, Any]) -> pd.DataFrame:
    df = read_table(source["path"], source.get("sheet_name"))
    if df.empty:
        return pd.DataFrame()

    triplet_size = int(source.get("triplet_size", 3))
    mean_offset = int(source.get("mean_offset", 0))
    dose_col = source.get("dose_column") or df.columns[0]
    ph = to_float(source.get("ph"))
    if ph is None:
        raise ValueError(
            f"Source {source['path']} uses triplet_excel format and needs a numeric 'ph' value."
        )

    extra_columns = list(df.columns[1:])
    if len(extra_columns) % triplet_size != 0:
        raise ValueError(
            f"{source['path']}: expected 1 dose column plus groups of {triplet_size} "
            f"columns (mean, SD, N). Got {len(extra_columns)} data columns."
        )

    records: list[dict[str, Any]] = []
    n_samples = len(extra_columns) // triplet_size

    for sample_index in range(n_samples):
        mean_col = extra_columns[sample_index * triplet_size + mean_offset]
        label = clean_excel_duplicate_suffix(mean_col)
        conc = source.get("conc")
        conc = to_float(conc) if conc is not None else parse_concentration_from_label(label)
        if conc is None:
            raise ValueError(
                f"Cannot find concentration in column '{mean_col}'. "
                "Use labels like 'Sample (0.05)' or set source['conc']."
            )

        material = clean_label(re.sub(r"\([^)]*\)", "", label))
        family = infer_family(material, source, config)

        for _, row in df.iterrows():
            dose = to_float(row[dose_col])
            value = to_float(row[mean_col])
            if dose is None or value is None:
                continue
            record = {
                "material": material,
                "family": family,
                "pH": ph,
                "dose": dose,
                "conc": conc,
                "lnDEF": value,
            }
            add_optional_metadata(record, row, source)
            records.append(record)

    return pd.DataFrame(records)


def guess_format(source: dict[str, Any]) -> str:
    df = read_table(source["path"], source.get("sheet_name"))
    lower_columns = {compact_key(col) for col in df.columns}
    has_long = (
        any(compact_key(alias) in lower_columns for alias in COLUMN_ALIASES["dose"])
        and any(compact_key(alias) in lower_columns for alias in COLUMN_ALIASES["conc"])
        and any(compact_key(alias) in lower_columns for alias in COLUMN_ALIASES["lnDEF"])
    )
    if has_long:
        return "long"

    if any(extract_dose_from_column(col) is not None for col in df.columns):
        return "wide_dose"

    if df.shape[1] > 1 and (df.shape[1] - 1) % int(source.get("triplet_size", 3)) == 0:
        return "triplet_excel"

    raise ValueError(
        f"Cannot detect data format for {source['path']}. "
        "Set source['format'] to 'long', 'wide_dose' or 'triplet_excel'."
    )


def load_source(source: dict[str, Any], config: dict[str, Any]) -> pd.DataFrame:
    path = Path(source["path"])
    if not path.exists():
        print(f"Missing source file, skipping: {path}")
        return pd.DataFrame()

    source_format = source.get("format", "auto")
    if source_format == "auto":
        source_format = guess_format(source)

    loaders = {
        "long": load_long_table,
        "wide": load_wide_dose_table,
        "wide_dose": load_wide_dose_table,
        "triplet": load_triplet_excel,
        "triplet_excel": load_triplet_excel,
    }
    if source_format not in loaders:
        raise ValueError(
            f"Unsupported format '{source_format}'. Use long, wide_dose or triplet_excel."
        )

    loaded = loaders[source_format](source, config)
    if not loaded.empty:
        include_families = source.get("include_families")
        exclude_families = source.get("exclude_families")
        if include_families:
            include = {clean_label(family) for family in include_families}
            loaded = loaded[loaded["family"].isin(include)]
        if exclude_families:
            exclude = {clean_label(family) for family in exclude_families}
            loaded = loaded[~loaded["family"].isin(exclude)]

    if loaded.empty:
        print(f"No usable rows loaded from: {path}")
    else:
        loaded["source"] = source.get("name", path.name)
    return loaded


def load_dataset(config: dict[str, Any] | None = None) -> pd.DataFrame:
    config = config or DEFAULT_CONFIG
    frames = [load_source(source, config) for source in config.get("sources", [])]
    frames = [frame for frame in frames if not frame.empty]
    if not frames:
        raise ValueError("No input data was loaded. Check file paths and source formats.")

    df = pd.concat(frames, ignore_index=True)
    group_cols = [col for col in ["material", "family", "pH", "dose", "conc"] if col in df.columns]
    for optional_col in ["experiment", "cell_line", "condition"]:
        if optional_col in df.columns:
            group_cols.append(optional_col)

    df = (
        df.groupby(group_cols, as_index=False)["lnDEF"]
        .mean()
        .sort_values(["family", "material", "pH", "dose", "conc"], kind="stable")
    )
    return df


def ordered_values(values: pd.Series, explicit_order: list[Any] | None = None) -> list[Any]:
    unique_values = list(pd.unique(values.dropna()))
    if not explicit_order:
        return sorted(unique_values, key=natural_key)

    order_as_text = [str(item) for item in explicit_order]
    by_text = {str(value): value for value in unique_values}
    ordered = [by_text[item] for item in order_as_text if item in by_text]
    ordered.extend(sorted([value for value in unique_values if str(value) not in order_as_text], key=natural_key))
    return ordered


def build_norm(values: pd.Series, center: float | None, symmetric: bool) -> Normalize:
    finite = pd.to_numeric(values, errors="coerce").dropna()
    if finite.empty:
        return Normalize(vmin=-1, vmax=1)

    vmin = float(finite.min())
    vmax = float(finite.max())
    if vmin == vmax:
        pad = abs(vmin) * 0.05 or 1.0
        vmin -= pad
        vmax += pad

    if center is None:
        return Normalize(vmin=vmin, vmax=vmax)

    if symmetric:
        limit = max(abs(vmin - center), abs(vmax - center))
        vmin = center - limit
        vmax = center + limit

    if not (vmin < center < vmax):
        return Normalize(vmin=vmin, vmax=vmax)
    return TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax)


def output_paths_for_group(group: str, plot_config: dict[str, Any]) -> list[Path]:
    output_dir = Path(plot_config.get("output_dir", "heatmap_outputs"))
    output_dir.mkdir(parents=True, exist_ok=True)
    output_names = plot_config.get("output_names", {})
    configured = output_names.get(group)

    if configured is None:
        return [output_dir / f"heatmap_{slugify(group)}.png"]
    if isinstance(configured, str):
        configured = [configured]
    return [output_dir / name for name in configured]


def plot_heatmaps(df: pd.DataFrame, plot_config: dict[str, Any] | None = None) -> list[Path]:
    plot_config = plot_config or {}
    if df.empty:
        raise ValueError("Cannot plot an empty dataframe.")

    group_column = plot_config.get("group_column", "family")
    row_column = plot_config.get("row_column", "material")
    panel_column = plot_config.get("panel_column", "pH")
    x_column = plot_config.get("x_column", "dose")
    y_column = plot_config.get("y_column", "conc")
    value_column = plot_config.get("value_column", "lnDEF")

    group_values = ordered_values(df[group_column]) if group_column in df.columns else ["All samples"]
    output_paths: list[Path] = []
    center = plot_config.get("center", 0)
    center = None if center in ("", None, "none", "None") else float(center)
    symmetric = bool(plot_config.get("symmetric_color_scale", True))
    global_scale = bool(plot_config.get("global_color_scale", True))
    annotate = bool(plot_config.get("annotate", True))
    cmap = plt.get_cmap(plot_config.get("cmap", "bwr")).copy()
    cmap.set_bad(plot_config.get("missing_color", "#f0f0f0"))
    global_norm = build_norm(df[value_column], center, symmetric) if global_scale else None

    for group in group_values:
        group_df = df[df[group_column] == group].copy() if group_column in df.columns else df.copy()
        if group_df.empty:
            continue

        material_order = plot_config.get("material_order", {}).get(str(group))
        rows = ordered_values(group_df[row_column], material_order)
        panels = ordered_values(group_df[panel_column])
        x_values = ordered_values(df[x_column] if global_scale else group_df[x_column])
        norm = global_norm or build_norm(group_df[value_column], center, symmetric)

        fig_height = max(2.6, float(plot_config.get("row_height", 2.7)) * len(rows))
        fig_width = max(6.0, float(plot_config.get("panel_width", 4.2)) * len(panels) + 1.4)
        fig, axes = plt.subplots(
            nrows=len(rows),
            ncols=len(panels),
            figsize=(fig_width, fig_height),
            squeeze=False,
            sharex=False,
        )

        last_image = None
        for row_index, material in enumerate(rows):
            for panel_index, panel in enumerate(panels):
                ax = axes[row_index, panel_index]
                subset = group_df[
                    (group_df[row_column] == material) & (group_df[panel_column] == panel)
                ]
                if subset.empty:
                    ax.axis("off")
                    continue

                pivot = (
                    subset.pivot_table(
                        index=y_column,
                        columns=x_column,
                        values=value_column,
                        aggfunc="mean",
                    )
                    .sort_index(axis=0)
                    .reindex(x_values, axis=1)
                )
                y_values = list(pivot.index)
                values = np.ma.masked_invalid(pivot.to_numpy(dtype=float))

                image = ax.imshow(values, cmap=cmap, norm=norm, aspect="auto")
                last_image = image

                ax.set_xticks(range(len(x_values)))
                ax.set_xticklabels([f"{value:g}" if isinstance(value, float) else value for value in x_values])
                if row_index == len(rows) - 1:
                    ax.set_xlabel(plot_config.get("dose_label", "Dose"))

                ax.set_yticks(range(len(y_values)))
                ax.set_yticklabels([f"{value:g}" if isinstance(value, float) else value for value in y_values])
                if panel_index == 0:
                    ax.set_ylabel(plot_config.get("concentration_label", "Concentration"))
                    ax.text(
                        -0.35,
                        0.5,
                        format_material_label(str(material)),
                        transform=ax.transAxes,
                        rotation=90,
                        va="center",
                        ha="center",
                        fontsize=12,
                        fontweight="bold",
                    )
                else:
                    ax.set_ylabel("")

                if annotate:
                    for y_pos, _ in enumerate(y_values):
                        for x_pos, _ in enumerate(x_values):
                            value = pivot.iat[y_pos, x_pos]
                            if pd.notna(value):
                                ax.text(x_pos, y_pos, f"{value:.2f}", ha="center", va="center", fontsize=8)

                if row_index == 0:
                    panel_label = plot_config.get("panel_label", panel_column)
                    ax.set_title(f"{panel_label} {panel}", fontsize=11)

        title = plot_config.get("group_titles", {}).get(str(group), str(group))
        fig.text(0.52, 0.965, title, ha="center", va="top", fontsize=14, fontweight="bold")
        fig.subplots_adjust(
            left=float(plot_config.get("left_margin", 0.22)),
            right=float(plot_config.get("right_margin", 0.82)),
            top=float(plot_config.get("top_margin", 0.92)),
            bottom=float(plot_config.get("bottom_margin", 0.08)),
            wspace=float(plot_config.get("wspace", 0.25)),
            hspace=float(plot_config.get("hspace", 0.55)),
        )

        if last_image is not None:
            cax = fig.add_axes([
                float(plot_config.get("colorbar_left", 0.86)),
                float(plot_config.get("colorbar_bottom", 0.15)),
                float(plot_config.get("colorbar_width", 0.03)),
                float(plot_config.get("colorbar_height", 0.70)),
            ])
            colorbar = fig.colorbar(last_image, cax=cax)
            label = plot_config.get("value_label", value_column)
            note = plot_config.get("colorbar_note", "")
            colorbar.set_label(f"{label}\n{note}" if note else label, fontsize=9)

        for output_path in output_paths_for_group(str(group), plot_config):
            fig.savefig(output_path, dpi=int(plot_config.get("dpi", 300)), bbox_inches="tight")
            print(f"Saved {group} heatmap to: {output_path}")
            output_paths.append(output_path)
        plt.close(fig)

    return output_paths


def load_config(path: str | Path | None) -> dict[str, Any]:
    if path is None:
        return DEFAULT_CONFIG
    with Path(path).open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    merged = DEFAULT_CONFIG.copy()
    merged.update(config)
    if "plot" in config:
        merged["plot"] = DEFAULT_CONFIG["plot"].copy()
        merged["plot"].update(config["plot"])
    return merged


def write_example_config(path: str | Path) -> None:
    path = Path(path)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(DEFAULT_CONFIG, handle, indent=2, ensure_ascii=False)
    print(f"Wrote example config to: {path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build adaptable lnDEF heatmaps from long CSV, wide dose CSV, "
            "or Excel triplet tables."
        )
    )
    parser.add_argument("--config", help="Path to a JSON config. Defaults to the bundled example setup.")
    parser.add_argument("--write-example-config", metavar="PATH", help="Write a reusable JSON config and exit.")
    parser.add_argument("--preview", action="store_true", help="Print the normalized dataset before plotting.")
    parser.add_argument("--output-dir", help="Override plot output directory.")
    parser.add_argument("--show", action="store_true", help="Show figures interactively after saving.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.write_example_config:
        write_example_config(args.write_example_config)
        return

    config = load_config(args.config)
    if args.output_dir:
        config.setdefault("plot", {})["output_dir"] = args.output_dir

    df = load_dataset(config)
    if args.preview:
        print(df.to_string(index=False))

    output_paths = plot_heatmaps(df, config.get("plot", {}))
    if args.show:
        for output_path in output_paths:
            print(output_path)
        plt.show()


if __name__ == "__main__":
    main()

from datetime import date
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# USER SETTINGS
# ============================================================
today = date.today()


filename = "/home/daniel/LibraFiles/CleanThesis/RootOutputs/4_20_2026_thesis_MB_von_mises_widths.txt"

# Choose observable:
#   "width"
#   "yield"
#   "baseline_B"
#   "ue_yield_counts"
#   "jet_yield_counts"
#   "ue_yield_fraction"
#   "jet_yield_fraction"
observable = "jet_yield_fraction"

# Choose peak for peak-specific observables: None, "near", or "away".
# Ignored for global observables such as the pedestal, UE yield, and jet yield.
peak_to_plot = "near"

# Choose scan direction:
#   "d0"       -> x-axis is D0 pT, anti-D0 pT is held constant
#   "anti_d0"  -> x-axis is anti-D0 pT, D0 pT is held constant
scan_axis = "d0"

# Save figures?
save_figures = True
output_dir = "/home/daniel/LibraFiles/CleanThesis/PostAnalysis"

# ------------------------------------------------------------
# Error-bar / quality control options
# ------------------------------------------------------------

# Remove points with NaN or nonpositive errors
drop_bad_errors = True

# Hard upper cut on absolute error (set to None to disable)
# For widths, units are rad. For yields, units are counts.
max_abs_error = None
# Example:
# max_abs_error = 0.8   # for widths
# max_abs_error = 2000  # for yields

# Relative error cut: remove points with err/value > this
# Set to None to disable
max_rel_error = None   # e.g. 1.0 = 100%

# If value is very close to zero, avoid division issues
small_value_epsilon = 1e-12

# ------------------------------------------------------------
# Plot range control
# ------------------------------------------------------------

# Option 1: determine y-range from central values only
use_central_values_for_ylim = True

# Add a margin around automatic y-range
ylim_padding_fraction = 0.15

# Optional manual y-limits
manual_ylim = None
# Example:
# manual_ylim = (0.0, 1.5)

# If you want percentile-based y-limits instead of min/max:
use_percentile_ylim = False
ylim_percentiles = (5, 95)

# ============================================================
# BIN DEFINITIONS
# ============================================================

PT_BINS = [
    ("first_pT", 0.0, 0.5),
    ("second_pT", 0.5, 1.0),
    ("third_pT", 1.0, 1.5),
    ("fourth_pT", 1.5, 2.0),
    ("fifth_pT", 2.0, 2.5),
    ("sixth_pT", 2.5, 3.0),
    ("seventh_pT", 3.0, 3.5),
    ("eighth_pT", 3.5, 8.0),
]
# PT_BINS = [
#     ("first_pT", 0.0, 1.0),
#     ("second_pT", 1.0, 2.0),
#     ("third_pT", 2.0, 3.0),
#     ("fourth_pT", 3.0, 4.0),
#     ("fifth_pT", 4.0, 5.0),
#     ("sixth_pT", 5.0, 6.0),
#     ("seventh_pT", 6.0, 7.0),
#     ("eighth_pT", 7.0, 8.0),
# ]


MULT_BINS = [
    ("low_mult", None, None),
    ("mid_mult", None, None),
    ("high_mult", None, None),
]

# ============================================================
# ORDERING / LABELS
# ============================================================

pt_order = [x[0] for x in PT_BINS]
pt_midpoints = {label: 0.5 * (low + high) for label, low, high in PT_BINS}
pt_ranges = {label: (low, high) for label, low, high in PT_BINS}
mult_order = [x[0] for x in MULT_BINS]

def pt_pretty_label(pt_class):
    low, high = pt_ranges[pt_class]
    # return f"{low:.0f}–{high:.0f} GeV/c"
    return f"{low:.1f}–{high:.1f} GeV/c"

def pt_axis_ticklabels():
    return [pt_pretty_label(label) for label in pt_order]

# def observable_columns(obs):
#     if obs == "width":
#         return "width_rad", "width_err_rad", "Width [rad]"
#     elif obs == "yield":
#         return "yield_counts", "yield_err_counts", "Yield [counts]"
#     else:
#         raise ValueError("observable must be 'width' or 'yield'")
def observable_columns(obs):
    if obs == "width":
        return "width_rad", "width_err_rad", "Width [rad]"
    elif obs == "yield":
        return "yield_counts", "yield_err_counts", "Yield [counts]"
    elif obs == "baseline_B":
        return "baseline_B", "baseline_B_err", "Pedestal height [counts/bin]"
    elif obs == "ue_yield_counts":
        return "ue_yield_counts", "ue_yield_err_counts", "UE Yield [counts]"
    elif obs == "jet_yield_counts":
        return "jet_yield_counts", "jet_yield_err_counts", "Jet Yield [counts]"
    elif obs == "ue_yield_fraction":
        return "ue_yield_fraction", "ue_yield_fraction_err", "UE Yield / histogram counts"
    elif obs == "jet_yield_fraction":
        return "jet_yield_fraction", "jet_yield_fraction_err", "Jet Yield / histogram counts"
    else:
        raise ValueError(
            "observable must be one of: "
            "'width', 'yield', 'baseline_B', 'ue_yield_counts', 'jet_yield_counts', "
            "'ue_yield_fraction', or 'jet_yield_fraction'"
        )

def scan_config(scan_axis):
    if scan_axis == "d0":
        return {
            "x_class_col": "d0_pt_class",
            "fixed_class_col": "anti_d0_pt_class",
            "x_mid_col_name": "d0_pt_mid",
            "x_label": r"$p_T(D^0)$ [GeV/$c$]",
            "fixed_label_tex": r"$p_T(\bar{D}^0)$",
            "filename_tag": "vs_d0pt_fixedAnti"
        }
    elif scan_axis == "anti_d0":
        return {
            "x_class_col": "anti_d0_pt_class",
            "fixed_class_col": "d0_pt_class",
            "x_mid_col_name": "anti_d0_pt_mid",
            "x_label": r"$p_T(\bar{D}^0)$ [GeV/$c$]",
            "fixed_label_tex": r"$p_T(D^0)$",
            "filename_tag": "vs_antid0pt_fixedD0"
        }
    else:
        raise ValueError("scan_axis must be 'd0' or 'anti_d0'")

# ============================================================
# READ FILE
# ============================================================

# If your file is actually whitespace-separated, change sep to r"\s+"
df = pd.read_csv(filename, sep="\t")

# Enforce categorical ordering
df["d0_pt_class"] = pd.Categorical(df["d0_pt_class"], categories=pt_order, ordered=True)
df["anti_d0_pt_class"] = pd.Categorical(df["anti_d0_pt_class"], categories=pt_order, ordered=True)
df["multiplicity_class"] = pd.Categorical(df["multiplicity_class"], categories=mult_order, ordered=True)

# Add midpoint columns
df["d0_pt_mid"] = df["d0_pt_class"].map(pt_midpoints).astype(float)
df["anti_d0_pt_mid"] = df["anti_d0_pt_class"].map(pt_midpoints).astype(float)

if "total_hist_counts" not in df.columns and "total_model_yield_counts" in df.columns:
    # Backward-compatible fallback for older output files. Regenerate the
    # analysis output to normalize by the true histogram count integral.
    df["total_hist_counts"] = df["total_model_yield_counts"]

if "total_hist_counts" in df.columns:
    denom = df["total_hist_counts"].replace(0, np.nan)
    if "ue_yield_fraction" not in df.columns and "ue_yield_counts" in df.columns:
        df["ue_yield_fraction"] = df["ue_yield_counts"] / denom
        df["ue_yield_fraction_err"] = df["ue_yield_err_counts"] / denom
    if "jet_yield_fraction" not in df.columns and "jet_yield_counts" in df.columns:
        df["jet_yield_fraction"] = df["jet_yield_counts"] / denom
        df["jet_yield_fraction_err"] = df["jet_yield_err_counts"] / denom

# ============================================================
# CHOOSE OBSERVABLE AND SCAN CONFIG
# ============================================================

value_col, err_col, y_label = observable_columns(observable)
cfg = scan_config(scan_axis)
global_observable = observable in {
    "baseline_B",
    "ue_yield_counts",
    "jet_yield_counts",
    "ue_yield_fraction",
    "jet_yield_fraction",
}

if global_observable:
    # These observables are written on both the near and away rows for the same
    # fit. Keep one copy so the plotter does not combine duplicate points.
    df = df[df["peak"] == "near"].copy()
elif peak_to_plot is not None:
    df = df[df["peak"] == peak_to_plot].copy()

x_class_col = cfg["x_class_col"]
fixed_class_col = cfg["fixed_class_col"]
x_mid_col_name = cfg["x_mid_col_name"]
x_label = cfg["x_label"]
fixed_label_tex = cfg["fixed_label_tex"]
filename_tag = cfg["filename_tag"]

# ============================================================
# FILTERING FUNCTION
# ============================================================

def filter_bad_points(subdf, value_col, err_col):
    out = subdf.copy()

    # Remove missing values
    out = out[out[value_col].notna() & out[err_col].notna()]

    if drop_bad_errors:
        out = out[np.isfinite(out[value_col]) & np.isfinite(out[err_col])]
        out = out[out[err_col] > 0]

    if max_abs_error is not None:
        out = out[out[err_col] <= max_abs_error]

    if max_rel_error is not None:
        denom = np.maximum(np.abs(out[value_col].to_numpy()), small_value_epsilon)
        rel_err = np.abs(out[err_col].to_numpy()) / denom
        out = out[rel_err <= max_rel_error]

    return out

# ============================================================
# DUPLICATE HANDLING
# ============================================================

def combine_duplicates(this_df, x_class_col, value_col, err_col):
    """
    If multiple rows appear for the same x-bin, combine them with inverse-variance weighting.
    """
    grouped_rows = []

    for xclass in pt_order:
        xbin = this_df[this_df[x_class_col] == xclass]
        if xbin.empty:
            continue

        if len(xbin) == 1:
            row = xbin.iloc[0]
            grouped_rows.append({
                x_class_col: xclass,
                x_mid_col_name: pt_midpoints[xclass],
                value_col: row[value_col],
                err_col: row[err_col],
            })
        else:
            vals = xbin[value_col].to_numpy(dtype=float)
            errs = xbin[err_col].to_numpy(dtype=float)

            good = errs > 0
            if np.any(good):
                weights = 1.0 / errs[good]**2
                avg = np.sum(weights * vals[good]) / np.sum(weights)
                err = np.sqrt(1.0 / np.sum(weights))
            else:
                avg = np.mean(vals)
                err = np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0.0

            grouped_rows.append({
                x_class_col: xclass,
                x_mid_col_name: pt_midpoints[xclass],
                value_col: avg,
                err_col: err,
            })

    return pd.DataFrame(grouped_rows)

# ============================================================
# Y-LIMIT HELPER
# ============================================================

def compute_ylims(all_y, all_yerr):
    if manual_ylim is not None:
        return manual_ylim

    if len(all_y) == 0:
        return None

    y = np.array(all_y, dtype=float)
    ye = np.array(all_yerr, dtype=float)

    if use_central_values_for_ylim:
        vals = y
    else:
        vals = np.concatenate([y - ye, y + ye])

    if use_percentile_ylim and len(vals) > 2:
        ymin, ymax = np.percentile(vals, ylim_percentiles)
    else:
        ymin, ymax = np.min(vals), np.max(vals)

    if np.isclose(ymin, ymax):
        pad = 0.1 * (abs(ymax) if ymax != 0 else 1.0)
        ymin -= pad
        ymax += pad
    else:
        pad = ylim_padding_fraction * (ymax - ymin)
        ymin -= pad
        ymax += pad

    return ymin, ymax

# ============================================================
# PLOTTING
# ============================================================

marker_map = {
    "low_mult": "o",
    "mid_mult": "s",
    "high_mult": "^",
}

linestyle_map = {
    "low_mult": "-",
    "mid_mult": "--",
    "high_mult": "-.",
}

if save_figures:
    os.makedirs(output_dir, exist_ok=True)

peaks_to_loop = ["global"] if global_observable else sorted(df["peak"].dropna().unique())

for peak in peaks_to_loop:
    peak_df = df.copy() if global_observable else df[df["peak"] == peak].copy()

    for fixed_pt_class in pt_order:
        sub = peak_df[peak_df[fixed_class_col] == fixed_pt_class].copy()
        sub = filter_bad_points(sub, value_col, err_col)

        if sub.empty:
            continue

        plt.figure(figsize=(9, 6), dpi=150)

        anything_plotted = False
        all_y = []
        all_yerr = []

        for mult_class in mult_order:
            this = sub[sub["multiplicity_class"] == mult_class].copy()
            this = this.sort_values(x_class_col)

            if this.empty:
                continue

            gdf = combine_duplicates(this, x_class_col, value_col, err_col)
            if gdf.empty:
                continue

            gdf = gdf.sort_values(x_mid_col_name)

            plt.errorbar(
                gdf[x_mid_col_name],
                gdf[value_col],
                yerr=gdf[err_col],
                fmt=marker_map[mult_class],
                linestyle=linestyle_map[mult_class],
                capsize=4,
                markersize=6,
                linewidth=1.5,
                label=mult_class.replace("_", " "),
            )

            all_y.extend(gdf[value_col].tolist())
            all_yerr.extend(gdf[err_col].tolist())
            anything_plotted = True

        if not anything_plotted:
            plt.close()
            continue

        fixed_label = pt_pretty_label(fixed_pt_class)

        plt.xlabel(x_label, fontsize=12)
        plt.ylabel(y_label, fontsize=12)
        if global_observable:
            title = f"{observable} vs pT\nFixed {fixed_label_tex}: {fixed_label}"
            output_peak_tag = "global"
        else:
            title = (
                f"{peak.capitalize()}-side {observable} vs pT\n"
                f"Fixed {fixed_label_tex}: {fixed_label}"
            )
            output_peak_tag = peak

        plt.title(title, fontsize=13)

        plt.xticks(
            [pt_midpoints[label] for label in pt_order],
            pt_axis_ticklabels(),
            rotation=45
        )

        ylims = compute_ylims(all_y, all_yerr)
        if ylims is not None:
            plt.ylim(*ylims)

        plt.legend(title="Multiplicity")
        plt.grid(alpha=0.3)
        plt.tight_layout()

        if save_figures:
            outname = (
                f"{output_dir}/"
                f"{today.month}_{today.day}_{today.year}_{observable}_{output_peak_tag}_{filename_tag}_{fixed_pt_class}.png"
            )
            plt.savefig(outname, bbox_inches="tight")

        plt.show()

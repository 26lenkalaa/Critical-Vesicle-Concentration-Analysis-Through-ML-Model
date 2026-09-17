"""
CVC (Critical Vesicle/Micelle Concentration) detection via supervised learning.

Trains a logistic regression classifier on labeled wells (above/below CVC)
from the two confirmed dilution series in this project, then applies it
back to those series as a sanity check.
"""

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

# ----------------------------- CONFIG -------------------------------
LOW_RANGE = (400, 530)
HIGH_RANGE = (400, 600)
BASELINE_WAVELENGTH = 600
OUTPUT_PLOT_PATH = "cvc_supervised.png"
# ----------------------------------------------------------------------


def load_spectrum(csv_path):
    """Load an absorbance spectrum CSV, skipping the instrument metadata header."""
    with open(csv_path, "r") as f:
        lines = f.readlines()
    header_row = next(
        i for i, line in enumerate(lines) if line.strip().startswith("Wavelength,")
    )
    df = pd.read_csv(csv_path, skiprows=header_row, low_memory=False)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(subset=["Wavelength"])
    df = df.dropna(axis=0, how="any")
    return df.reset_index(drop=True)


def find_sample_peaks(df, sample_col, low_range=LOW_RANGE, high_range=HIGH_RANGE):
    wl = df["Wavelength"]
    low_mask = (wl >= low_range[0]) & (wl <= low_range[1])
    high_mask = (wl >= high_range[0]) & (wl <= high_range[1])
    low_idx = df.loc[low_mask, sample_col].idxmax()
    high_idx = df.loc[high_mask, sample_col].idxmax()
    return {"low_value": df.loc[low_idx, sample_col], "high_value": df.loc[high_idx, sample_col]}


def baseline_value(df, sample_col, baseline_wavelength=BASELINE_WAVELENGTH):
    row = df.loc[df["Wavelength"] == baseline_wavelength]
    if row.empty:
        row = df.iloc[[-1]]
    return row[sample_col].iloc[0]


def compute_well_ratios(df):
    sample_cols = [c for c in df.columns if c != "Wavelength"]
    ratios = {}
    for col in sample_cols:
        peaks = find_sample_peaks(df, col)
        baseline = baseline_value(df, col)
        high_sub = peaks["high_value"] - baseline
        low_sub = peaks["low_value"] - baseline
        ratios[col] = high_sub / low_sub if low_sub != 0 else float("nan")
    return ratios


# --------------------------------------------------------------------
# Labeled training data
#
# Wells listed lowest to highest concentration (index position stands in
# for concentration, since the real top standard concentration and
# spacing were never confirmed -- see the project README). label = 0
# below the confirmed breakpoint, 1 at or above it.
# --------------------------------------------------------------------

TRAINING_SERIES = [
    {
        "csv": "Sample26_Absorbance_Spectrum.csv",
        "wells_low_to_high": ["Un0014 (B05)", "Un0013 (B04)", "Un0012 (B03)", "Un0011 (B02)",
                               "Un0010 (B01)", "Un0021 (A12)", "Un0020 (A11)", "Un0019 (A10)",
                               "Un0009 (A09)", "Un0008 (A08)"],
        "breakpoint_index": 5,  # confirmed: between index 4 (Un0010) and index 5 (Un0021)
        "confidence": "high",
    },
    {
        "csv": "Sample15_Absorbance_Spectrum.csv",
        "wells_low_to_high": ["Un0013 (H01)", "Un0012 (G12)", "Un0011 (G11)", "Un0010 (G10)",
                               "Un0009 (G09)", "Un0008 (G08)", "Un0007 (G07)", "Un0006 (G06)",
                               "Un0005 (G05)", "Un0004 (G04)"],
        "breakpoint_index": 2,  # confirmed: between index 1 (Un0012) and index 2 (Un0011)
        "confidence": "low",   # noisy at the bottom -- see project README, Open Items
    },
]


def build_training_set():
    """Returns (X, y, meta) where X is [ratio, index_position] per well,
    y is the above/below CVC label, and meta records provenance for
    anyone auditing where a label came from."""
    X, y, meta = [], [], []
    for series in TRAINING_SERIES:
        df = load_spectrum(series["csv"])
        ratios = compute_well_ratios(df)
        for i, well in enumerate(series["wells_low_to_high"]):
            label = 0 if i < series["breakpoint_index"] else 1
            X.append([ratios[well], i])
            y.append(label)
            meta.append({"csv": series["csv"], "well": well, "index": i, "label": label,
                         "source_confidence": series["confidence"]})
    return np.array(X, dtype=float), np.array(y, dtype=int), meta


def train_classifier(X, y):
    """Logistic regression on standardized features. With only 20 labeled
    points, leave-one-out cross-validation is used to get an honest,
    if noisy, estimate of how well this generalizes -- a single train/test
    split would be too small to mean much either way."""
    model = make_pipeline(StandardScaler(), LogisticRegression())
    loo_scores = cross_val_score(model, X, y, cv=LeaveOneOut())
    model.fit(X, y)
    return model, loo_scores


def predict_breakpoint(model, ratios_low_to_high):
    """Given a new series' ratios (ordered lowest to highest concentration),
    predict per-well labels and report where the model's prediction flips
    from below-CVC to above-CVC."""
    X_new = np.array([[r, i] for i, r in enumerate(ratios_low_to_high)], dtype=float)
    probs = model.predict_proba(X_new)[:, 1]
    preds = (probs >= 0.5).astype(int)
    breakpoint_idx = None
    for i in range(len(preds) - 1):
        if preds[i] == 0 and preds[i + 1] == 1:
            breakpoint_idx = i + 1
            break
    return preds, probs, breakpoint_idx


def plot_training_data(X, y, meta, model, out_path=OUTPUT_PLOT_PATH):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot.")
        return

    markers = {"Sample26_Absorbance_Spectrum.csv": "o", "Sample15_Absorbance_Spectrum.csv": "s"}
    colors = {0: "tab:blue", 1: "tab:red"}

    plt.figure(figsize=(7, 5))
    for m, (ratio, idx) in zip(meta, X):
        plt.scatter(idx, ratio, marker=markers[m["csv"]], color=colors[m["label"]],
                    s=70, edgecolor="black")

    idx_range = np.linspace(X[:, 1].min(), X[:, 1].max(), 200)
    ratio_range = np.linspace(X[:, 0].min(), X[:, 0].max(), 200)
    grid_idx, grid_ratio = np.meshgrid(idx_range, ratio_range)
    grid_X = np.column_stack([grid_ratio.ravel(), grid_idx.ravel()])
    grid_probs = model.predict_proba(grid_X)[:, 1].reshape(grid_idx.shape)
    plt.contour(grid_idx, grid_ratio, grid_probs, levels=[0.5], colors="black", linestyles="--")

    plt.xlabel("Well index (ascending concentration)")
    plt.ylabel("Absorbance ratio")
    plt.title("Supervised CVC classifier: training wells and decision boundary\n"
              "circle=Sample26 (high-confidence labels), square=Sample15 (low-confidence labels)\n"
              "blue=below CVC, red=above CVC")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"\nPlot saved to {out_path}")


def main():
    print("Building labeled training set from the two confirmed series...")
    X, y, meta = build_training_set()
    print(f"  {len(y)} labeled wells ({sum(y)} above CVC, {len(y) - sum(y)} below)")
    for m in meta:
        print(f"    {m['csv']:35s} {m['well']:15s} index={m['index']} label={m['label']} "
              f"(source confidence: {m['source_confidence']})")

    model, loo_scores = train_classifier(X, y)
    print(f"\nLeave-one-out CV accuracy: {loo_scores.mean():.2f} "
          f"({int(loo_scores.sum())}/{len(loo_scores)} correct)")

    print("\nPredicted breakpoint on each training series:")
    for series in TRAINING_SERIES:
        df = load_spectrum(series["csv"])
        ratios = compute_well_ratios(df)
        ratio_seq = [ratios[w] for w in series["wells_low_to_high"]]
        preds, probs, bp_idx = predict_breakpoint(model, ratio_seq)
        print(f"  {series['csv']}: predicted index {bp_idx} "
              f"(confirmed index {series['breakpoint_index']})")

    plot_training_data(X, y, meta, model)


if __name__ == "__main__":
    main()
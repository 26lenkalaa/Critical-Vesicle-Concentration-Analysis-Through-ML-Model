"""
End-to-end Critical Vesicle/Micelle Concentration (CVC) automation.

Pipeline: raw plate-reader CSV -> per-well absorbance ratio -> per-well
known concentration -> automated breakpoint (CVC) detection -> plot + value.

No manual peak-picking or curve-reading is required once the CONFIG
section below is filled in with your actual serial dilution parameters.

--------------------------------------------------------------------
WHAT STILL NEEDS TO BE CONFIRMED BEFORE THIS OUTPUT IS TRUSTWORTHY
--------------------------------------------------------------------
The script currently assumes:
  1. every non-blank column in the file is one point in a single serial
     dilution series (not separate replicates or unrelated samples)
  2. the first column in the file is the highest concentration and each
     subsequent column is diluted by DILUTION_FACTOR from the one before it
  3. Blank1 is true concentration 0 (buffer/dye only, no lipid/surfactant)
  4. TOP_CONCENTRATION_MM below is your actual top standard's concentration

If any of those don't match your protocol, fix the CONFIG section (or the
column-ordering logic in `assign_concentrations`) before trusting the
printed CVC value. Everything downstream of CONFIG is fully automated.
--------------------------------------------------------------------

Two independent methods are used to detect the breakpoint, and both are
printed so you can sanity-check that they roughly agree:

  1. K-means clustering (unsupervised ML) on the ratio values, splitting
     wells into a "below CVC" cluster (baseline ratio) and an "above CVC"
     cluster (shifted ratio). The CVC estimate is the midpoint, in
     concentration space, between the last low-cluster point and the
     first high-cluster point going down the dilution series.

  2. A four-parameter sigmoid (logistic) fit of ratio vs. log10(concentration),
     the standard approach in the CMC/CVC literature. The fitted midpoint
     parameter, converted back out of log space, is the second CVC estimate.
"""

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from scipy.optimize import curve_fit

# ----------------------------- CONFIG -------------------------------
CSV_PATH = "Sample15_Absorbance_Spectrum.csv"
TOP_CONCENTRATION_MM = 10.0   # TODO: replace with your real top standard concentration
DILUTION_FACTOR = 2.0         # TODO: confirm serial dilution factor (2.0 = two-fold)
BLANK_LABEL_CONTAINS = "Blank"  # column names containing this are treated as concentration 0
LOW_RANGE = (400, 530)        # "low peak" search window, nm
HIGH_RANGE = (400, 600)       # "high peak" search window, nm
BASELINE_WAVELENGTH = 600     # nm, used for baseline subtraction
OUTPUT_PLOT_PATH = "cvc_curve.png"
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
    df = df.reset_index(drop=True)
    return df


def find_sample_peaks(df, sample_col, low_range=LOW_RANGE, high_range=HIGH_RANGE):
    """Find one sample column's own low-range and high-range peak."""
    wl = df["Wavelength"]
    low_mask = (wl >= low_range[0]) & (wl <= low_range[1])
    high_mask = (wl >= high_range[0]) & (wl <= high_range[1])

    low_idx = df.loc[low_mask, sample_col].idxmax()
    high_idx = df.loc[high_mask, sample_col].idxmax()

    return {
        "low_value": df.loc[low_idx, sample_col],
        "high_value": df.loc[high_idx, sample_col],
    }


def baseline_value(df, sample_col, baseline_wavelength=BASELINE_WAVELENGTH):
    row = df.loc[df["Wavelength"] == baseline_wavelength]
    if row.empty:
        row = df.iloc[[-1]]
    return row[sample_col].iloc[0]


def compute_well_ratios(df):
    """Return {well_name: ratio} for every sample column in the file."""
    sample_cols = [c for c in df.columns if c != "Wavelength"]
    ratios = {}
    for col in sample_cols:
        peaks = find_sample_peaks(df, col)
        baseline = baseline_value(df, col)
        high_sub = peaks["high_value"] - baseline
        low_sub = peaks["low_value"] - baseline
        ratios[col] = high_sub / low_sub if low_sub != 0 else float("nan")
    return ratios


def assign_concentrations(well_names, top_concentration, dilution_factor, blank_contains):
    """Assign a concentration to each well under a single serial dilution series.

    First non-blank column = top_concentration; each subsequent non-blank
    column is diluted by dilution_factor from the previous one. Any column
    whose name contains blank_contains is assigned concentration 0.
    """
    concentrations = {}
    dilution_step = 0
    for name in well_names:
        if blank_contains in name:
            concentrations[name] = 0.0
        else:
            concentrations[name] = top_concentration / (dilution_factor ** dilution_step)
            dilution_step += 1
    return concentrations


def detect_cvc_kmeans(concentrations, ratios):
    """Unsupervised breakpoint detection via 2-cluster K-means on ratio values."""
    wells = list(concentrations.keys())
    conc = np.array([concentrations[w] for w in wells])
    ratio = np.array([ratios[w] for w in wells])

    order = np.argsort(-conc)  # descending concentration
    conc, ratio = conc[order], ratio[order]

    km = KMeans(n_clusters=2, n_init=10, random_state=0)
    labels = km.fit_predict(ratio.reshape(-1, 1))

    # the cluster with the higher mean ratio is "above CVC"
    high_cluster = np.argmax(km.cluster_centers_.ravel())
    is_high = labels == high_cluster

    if is_high.all() or (~is_high).all():
        return None  # no transition found; every well fell in one cluster

    # walk down the (descending) concentration series and find where it
    # switches from "high" to "low" cluster membership
    for i in range(len(conc) - 1):
        if is_high[i] and not is_high[i + 1]:
            return (conc[i] + conc[i + 1]) / 2.0
    return None


def _sigmoid(log_conc, top, bottom, midpoint_log, slope):
    return bottom + (top - bottom) / (1.0 + np.exp(-slope * (log_conc - midpoint_log)))


def detect_cvc_sigmoid(concentrations, ratios):
    """Breakpoint via a fitted 4-parameter sigmoid of ratio vs log10(concentration)."""
    wells = [w for w in concentrations if concentrations[w] > 0]  # log10 needs > 0
    conc = np.array([concentrations[w] for w in wells])
    ratio = np.array([ratios[w] for w in wells])
    log_conc = np.log10(conc)

    p0 = [ratio.max(), ratio.min(), np.median(log_conc), 1.0]
    try:
        popt, _ = curve_fit(_sigmoid, log_conc, ratio, p0=p0, maxfev=10000)
    except RuntimeError:
        return None

    midpoint_log = popt[2]
    return 10 ** midpoint_log


def main(csv_path=CSV_PATH):
    df = load_spectrum(csv_path)
    ratios = compute_well_ratios(df)
    well_names = list(ratios.keys())
    concentrations = assign_concentrations(
        well_names, TOP_CONCENTRATION_MM, DILUTION_FACTOR, BLANK_LABEL_CONTAINS
    )

    print(f"Loaded {len(well_names)} wells from {csv_path}\n")
    for name in well_names:
        print(f"  {name}: concentration = {concentrations[name]:.6g}, ratio = {ratios[name]:.4f}")

    cvc_kmeans = detect_cvc_kmeans(concentrations, ratios)
    cvc_sigmoid = detect_cvc_sigmoid(concentrations, ratios)

    print("\nAutomated CVC estimates:")
    print(f"  K-means clustering breakpoint: {cvc_kmeans}")
    print(f"  Sigmoid fit inflection point:  {cvc_sigmoid}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        conc_sorted = sorted(concentrations.values())
        conc_arr = np.array([concentrations[w] for w in well_names])
        ratio_arr = np.array([ratios[w] for w in well_names])

        plt.figure(figsize=(7, 5))
        plt.scatter(conc_arr, ratio_arr, label="wells")
        if cvc_kmeans:
            plt.axvline(cvc_kmeans, color="orange", linestyle="--", label=f"K-means CVC = {cvc_kmeans:.4g}")
        if cvc_sigmoid:
            plt.axvline(cvc_sigmoid, color="green", linestyle=":", label=f"Sigmoid CVC = {cvc_sigmoid:.4g}")
        plt.xscale("log")
        plt.xlabel("Concentration (log scale)")
        plt.ylabel("High/low absorbance ratio")
        plt.title("CVC determination")
        plt.legend()
        plt.tight_layout()
        plt.savefig(OUTPUT_PLOT_PATH, dpi=150)
        print(f"\nPlot saved to {OUTPUT_PLOT_PATH}")
    except ImportError:
        print("\nmatplotlib not available; skipping plot.")


if __name__ == "__main__":
    main()
"""
Absorbance spectrum analysis for Sample15_Absorbance_Spectrum.csv

For each sample column in the plate:
    - find that sample's own peak absorbance across the full recorded
      wavelength range ("high peak")
    - find that sample's own peak absorbance in the lower sub-range,
      roughly 400-530 nm ("low peak")
    - baseline-subtract each peak against the absorbance at the last
      recorded wavelength (600 nm)
    - compute a high/low ratio for that sample

The last three sample columns are treated as known standards (default
concentrations 0, 5, 10) and used to fit a calibration line of
concentration vs. ratio. That line is then used to back-calculate an
estimated concentration for every other sample on the plate.

Compared to the original script, this version:
    - actually skips the instrument metadata rows instead of relying on
      dropna() to filter them out by accident
    - derives the wavelength range and column count from the data
      itself instead of hardcoding 400/201/23
    - finds each sample's peak independently instead of only checking
      one column and applying that wavelength to every sample
    - is explicit about which wavelength range counts as "high" vs "low"
"""

import numpy as np
import pandas as pd


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


def get_concentration_info():
    """Ask for the original concentration and stock percentage, return mmol value.

    NOTE: this is currently informational only. The calibration curve below
    uses the fixed standard_concentrations you pass to main(), not this
    value. Wire this in if the top standard's concentration should be
    derived from the stock dilution instead of hardcoded.
    """
    og_conc = float(input("What was the original concentration of the substance in moles: "))
    og_percent = float(input("Enter percentage of that concentration in your stock solution (no %% symbol): "))
    conc_in_mmol = (og_percent / 100.0) * og_conc * 1000
    print(f"Stock concentration: {conc_in_mmol:.4f} mmol")
    return conc_in_mmol


def find_sample_peaks(df, sample_col, low_range=(400, 530), high_range=(400, 600)):
    """Find each sample's own low-range and high-range peak."""
    wl = df["Wavelength"]

    low_mask = (wl >= low_range[0]) & (wl <= low_range[1])
    high_mask = (wl >= high_range[0]) & (wl <= high_range[1])

    low_idx = df.loc[low_mask, sample_col].idxmax()
    high_idx = df.loc[high_mask, sample_col].idxmax()

    return {
        "low_wavelength": wl[low_idx],
        "low_value": df.loc[low_idx, sample_col],
        "high_wavelength": wl[high_idx],
        "high_value": df.loc[high_idx, sample_col],
    }


def baseline_value(df, sample_col, baseline_wavelength=600):
    """Absorbance at the baseline wavelength (falls back to the last row if missing)."""
    row = df.loc[df["Wavelength"] == baseline_wavelength]
    if row.empty:
        row = df.iloc[[-1]]
    return row[sample_col].iloc[0]


def analyze_samples(df):
    """Compute peaks, baseline subtraction, and ratio for every sample column."""
    sample_cols = [c for c in df.columns if c != "Wavelength"]
    results = {}

    for col in sample_cols:
        peaks = find_sample_peaks(df, col)
        baseline = baseline_value(df, col)

        high_sub = peaks["high_value"] - baseline
        low_sub = peaks["low_value"] - baseline
        ratio = high_sub / low_sub if low_sub != 0 else float("nan")

        results[col] = {**peaks, "baseline": baseline, "high_sub": high_sub,
                         "low_sub": low_sub, "ratio": ratio}

    return results


def fit_standard_curve(results, standard_cols, standard_concentrations):
    """Fit concentration vs. ratio using the given standard columns."""
    if len(standard_cols) != len(standard_concentrations):
        raise ValueError("Number of standard columns must match number of known concentrations.")

    y = np.array([results[c]["ratio"] for c in standard_cols], dtype=float)
    x = np.array(standard_concentrations, dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    return slope, intercept


def main(csv_path, standard_concentrations=(0, 5, 10), n_standards=3):
    df = load_spectrum(csv_path)
    sample_cols = [c for c in df.columns if c != "Wavelength"]
    standard_cols = sample_cols[-n_standards:]

    results = analyze_samples(df)
    slope, intercept = fit_standard_curve(results, standard_cols, standard_concentrations)

    print(f"Loaded {len(df)} wavelength points ({df['Wavelength'].min()}-{df['Wavelength'].max()} nm) "
          f"across {len(sample_cols)} sample columns.\n")

    for col, r in results.items():
        print(f"{col}: high peak {r['high_value']:.4f} at {r['high_wavelength']} nm, "
              f"low peak {r['low_value']:.4f} at {r['low_wavelength']} nm, "
              f"ratio = {r['ratio']:.4f}")

    print(f"\nCalibration curve using standards {standard_cols} "
          f"at concentrations {standard_concentrations}:")
    print(f"  slope = {slope:.6f}, intercept = {intercept:.6f}")

    print("\nEstimated concentrations for remaining samples:")
    for col in sample_cols:
        if col in standard_cols:
            continue
        est_conc = (results[col]["ratio"] - intercept) / slope
        print(f"  {col}: {est_conc:.4f}")


if __name__ == "__main__":
    csv_file = input("Enter file name: ").strip() or "Sample15_Absorbance_Spectrum.csv"
    main(csv_file)
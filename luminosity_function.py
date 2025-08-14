"""
Module for computing the cumulative luminosity-weighted GAMA luminosity function.
"""

import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad


def build_integrated_lf(path: str = "gama_lf/lf.dat", cut: float = -14) -> interp1d:
    """
    Reads the GAMA luminosity function file, filters by magnitude, builds the luminosity-weighted
    luminosity function, integrates it, and returns an interpolation function over that integral.

    Args:
        path: Path to the LF file.
        cut: Magnitude cut (default = -14)

    Returns:
        interp1d: Interpolation function over the integrated phi(L) from -30 to M.
    """
    df = pd.read_csv(path, delim_whitespace=True)
    df = df[df["mag_bin_centre"] < cut].copy()

    mags = df["mag_bin_centre"].values
    phi = df["luminosity_function"].values

    # Compute luminosity-weighted φ: phi * 10^(-0.4 * M)
    phi_lum = phi * 10 ** (-0.4 * mags)

    # Interpolation of φ(L)
    func_lum = interp1d(mags, phi_lum,bounds_error=False,fill_value=(phi_lum[0], phi_lum[-1]))

    # Integrate from -30 to each M
    integrals = np.array([
        quad(func_lum, -30, M, limit=1000)[0] if np.isfinite(M) else 0.0
        for M in mags
    ])

    # Replace zeros with min positive value
    positive = integrals[integrals > 0]
    if len(positive) == 0:
        raise ValueError("All integrated values are zero; check your LF input.")
    min_val = positive.min()
    integrals[integrals == 0] = min_val

    return interp1d(mags, integrals, bounds_error=False, fill_value=(integrals[0], integrals[-1]))


if __name__ == '__main__':
    func = build_integrated_lf()
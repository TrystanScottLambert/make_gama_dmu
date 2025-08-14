"""
Module testing the luminosity function is correct for GAMA.
"""

import unittest
import numpy.testing as npt

from luminosity_function import build_integrated_lf

class LFTest(unittest.TestCase):
    """Testing the LF by comparing to the R code."""

    def test_lf(self):
        """Copmaring to the LFswmlintfunc in the old R code"""
        mags = [-18, -19, -20, -21, -22, -25]
        ans = [2067784.45, 1786564.75, 1209171.74, 394986.29, 27748.48, 1023.67]
        func = build_integrated_lf()
        res = func(mags)
        npt.assert_allclose(res, ans, rtol=1e-4)


if __name__ == '__main__':
    unittest.main()

"""
Unit tests for the ke correction
"""

import unittest
import numpy.testing as npt
from main import calc_ke_correction

class KETest(unittest.TestCase):
    """
    Testing that the KE correction matches the old R code.
    """
    def test_comparison_to_r(self):
        """Ran the KECorr function to get the answer that we are comparing to."""

        redshifts = [0.01, 0.02, 0.03, 0.10, 0.50, 1.00]
        result = calc_ke_correction(redshifts)

        ans = [
            -0.005927812,
            -0.01205652,
            -0.01837549,
            -0.06689517,
            -0.1963624,
            2.426494,
        ]
        npt.assert_almost_equal(result, ans, decimal=4)

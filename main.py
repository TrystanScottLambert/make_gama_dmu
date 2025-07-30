"""
Main script for building the GAMA Group Catalog
"""

from dataclasses import dataclass
import pandas as pd
import numpy as np
from nessie import RedshiftCatalog, FlatCosmology
from nessie.helper_funcs import create_density_function
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u


cosmo = FlatCosmology(0.7, 0.3)
astropy_cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
OVERFACTOR = 400  # The factor the randoms were expanded.
B0, R0 = 0.06, 32


def convert_jansky_to_ab(janksy_array: np.ndarray[float]) -> np.ndarray[float]:
    """
    Converts the flux in Jansky into AB apparent magnitudes.
    """
    return 2.5 * (23 - np.log10(janksy_array)) - 48.6


def read_in_input_cats(
    catalog_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Reads in the GAMA catalog returns the four fields as their own data frames
    """
    input_df = pd.read_csv(catalog_name)
    g09_df_input = input_df[(input_df["RAcen"] < 150) & (input_df["Z"] < 0.5)]
    g12_df_input = input_df[
        (input_df["RAcen"] > 150) & (input_df["RAcen"] < 200) & (input_df["Z"] < 0.5)
    ]
    g15_df_input = input_df[
        (input_df["RAcen"] > 200) & (input_df["RAcen"] < 250) & (input_df["Z"] < 0.5)
    ]
    g23_df_input = input_df[(input_df["RAcen"] > 300) & (input_df["Z"] < 0.5)]
    return g09_df_input, g12_df_input, g15_df_input, g23_df_input


def calc_ke_correction(redshifts: np.ndarray[float]):
    """
    calculates the k+e correction for a given redshifts
    """
    redshift = np.asarray(redshifts)
    kcorrvals = np.array([0.20848, 1.0226, 0.52366, 3.5902, 2.3843])
    powers = (redshift[..., None] - 0.2) ** np.arange(len(kcorrvals))
    k = powers @ kcorrvals
    return k - 1.75 * redshift


@dataclass
class Field:
    """
    Stores the information for each gama field and generates the redshift catalog object in the
    post __init__
    """

    input_data_frame: pd.DataFrame
    fractional_area: float
    name: str

    def __post_init__(self):
        obs_df = self.input_data_frame[self.input_data_frame["NQ"] > 2]
        obs_df["apparent_mags"] = convert_jansky_to_ab(obs_df["flux_rl"])
        obs_df["absolute_mags"] = (
            obs_df["apparent_mags"]
            - cosmo.dist_mod(obs_df["Z"])
            - calc_ke_correction(obs_df["Z"])
        )

        randoms = np.loadtxt(f"randoms/gama_{self.name}_randoms.txt", skiprows=1)
        rho = create_density_function(
            randoms, len(randoms) / OVERFACTOR, self.fractional_area, cosmo
        )
        redcat = RedshiftCatalog(
            obs_df["RAcen"], obs_df["Deccen"], obs_df["Z"], rho, cosmo
        )

        search_radii = astropy_cosmo.arcsec_per_kpc_comoving(obs_df["Z"]) * 1 * u.Mpc
        search_radii_deg = search_radii.to(u.deg).value
        search_radii_deg[search_radii_deg > 180] = (
            180.0  # for values at z=0 essentially
        )

        redcat.calculate_completeness(
            self.input_data_frame["RAcen"],
            self.input_data_frame["Deccen"],
            search_radii_deg,
        )
        redcat.completeness = np.array(redcat.completeness)
        redcat.run_fof(B0, R0)

        self.redshift_catalog = redcat
        self.obs_df = obs_df


if __name__ == "__main__":
    INPUT_CAT = "gama_catalogs/gama_input_galaxies.csv"
    g09_input, g12_input, g15_input, g23_input = read_in_input_cats(INPUT_CAT)

    field_strings = ["g09", "g12", "g15", "g23"]
    fractional_areas = {
        "g09": 0.001453924,
        "g12": 0.001453924,
        "g15": 0.001453924,
        "g23": 0.001226274,
    }

    g09 = Field(g09_input, fractional_areas["g09"], "g09")
    g12 = Field(g12_input, fractional_areas["g12"], "g12")
    g15 = Field(g15_input, fractional_areas["g15"], "g15")
    g23 = Field(g23_input, fractional_areas["g23"], "g23")

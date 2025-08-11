"""
Main script for building the GAMA Group Catalog
"""

from dataclasses import dataclass
import pandas as pd
import numpy as np
from nessie import RedshiftCatalog
from nessie.helper_funcs import create_density_function
import astropy.units as u

from config import cosmo, astropy_cosmo, GROUP_ID_OFFSET, APPARENT_MAG_LIM, OVERFACTOR, B0, R0, SUN_MAG, MASS_A, MASS_FUNC_PARAMS, LUM_B, LUM_FUNC_PARAMS, AB_CUT
from luminosity_function import build_integrated_lf
from separations import add_separation_metrics

def functional_correction(multiplicity: np.ndarray[float], median_redshift: np.ndarray[float], params: np.ndarray[float]) -> np.ndarray[float]:
    """
    This is the general form of the plane equation used in equations 19 and 23 in R11.

    multiplicity of the group
    median_redshift of the group

    Returns the functional correction to be applied to the mass proxy.
    """
    return params[0] + params[1] * (multiplicity**(-0.5)) + params[2]*(median_redshift**(-0.5))

def convert_jansky_to_apparent(janksy_array: np.ndarray[float]) -> np.ndarray[float]:
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
        obs_df["apparent_mags"] = convert_jansky_to_apparent(obs_df["flux_rl"])
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

    def get_group_dmu(self) -> pd.DataFrame:
        """Creates the GAMA DMU for the groups."""
        vel_errors = np.repeat(50, len(self.obs_df))
        properties = pd.DataFrame(self.redshift_catalog.calculate_group_table(self.obs_df['absolute_mags'], vel_errors))
        group_ob_limit = APPARENT_MAG_LIM - cosmo.dist_mod(properties['median_redshift'])
        int_function = build_integrated_lf()
        lum_factor = int_function(AB_CUT)/int_function(group_ob_limit)
        properties['lum_corrected_mass'] = properties['mass_proxy'] * lum_factor
        properties['lum_corrected_flux'] = properties['flux_proxies'] * lum_factor
        properties['MassA'] = properties['lum_corrected_mass'] * MASS_A
        properties['LumB'] = properties['lum_corrected_flux'] * LUM_B * 10**(0.4 * SUN_MAG)
        properties['MassAfunc'] = properties['lum_corrected_mass'] * functional_correction(properties['multiplicity'], properties['median_redshift'], MASS_FUNC_PARAMS)
        properties['LumBfunc'] = properties['lum_corrected_flux'] * functional_correction(properties['multiplicity'], properties['median_redshift'], LUM_FUNC_PARAMS)
        properties['group_id'] = np.array(properties['group_id'] + GROUP_ID_OFFSET[self.name]).astype(int)
        properties['iter_uber_id'] = np.array(self.obs_df['UberID'])[np.array(properties['iter_idx'])]
        properties['bcg_uber_id'] = np.array(self.obs_df['UberID'])[np.array(properties['bcg_idxs'])]
        return properties

    def get_pair_dmu(self) -> pd.DataFrame:
        """Creates the GAMA DMU for the pairs."""
        pair_properties = pd.DataFrame(self.redshift_catalog.calculate_pair_table(self.obs_df['absolute_mags']))
        pair_properties['uber_id_1'] = np.array(self.obs_df['UberID'])[np.array(pair_properties['idx_1'])]
        pair_properties['uber_id_2'] = np.array(self.obs_df['UberID'])[np.array(pair_properties['idx_2'])]
        return pair_properties

    def get_galaxy_dmu(self) -> pd.DataFrame:
        """Creates the GAMA DMU for the galaxies."""
        new_dmu = self.obs_df.copy()
        new_dmu['group_ids'] = np.array(self.redshift_catalog.group_ids + GROUP_ID_OFFSET[self.name]).astype(int)
        new_dmu.loc[new_dmu["group_ids"] == 99999, "group_ids"] = 0 # zero for ungrouped.
        return new_dmu
 
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

    groups = g09.get_group_dmu()
    pairs = g09.get_pair_dmu()
    galaxies = g09.get_galaxy_dmu()

    print('Adding separation metrics.')
    galaxies = add_separation_metrics(galaxies, groups)

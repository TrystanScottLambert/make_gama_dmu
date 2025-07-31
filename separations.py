"""
Module to calcualte the separation metrics for the galaxy tables.
"""

from dataclasses import dataclass
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import pandas as pd
from config import cosmo


@dataclass
class Separations:
    """
    Class storing the positional information for a group
    """

    center: SkyCoord
    galaxy_positions: SkyCoord
    group_redshift: np.ndarray[float]
    name: str

    def __post_init__(self) -> None:
        """Calculate the angular separation and comoving distance. Used throughout"""
        self.ang_sep = self.center.separation(self.galaxy_positions).to(u.rad)
        self.distance = cosmo.comoving_distance(self.group_redshift)

    @property
    def angular_separation(self) -> np.ndarray[float]:
        """The angular separation between every galaxy and the center in arcseconds."""
        return np.array(self.ang_sep.to(u.arcsec).value)

    @property
    def projected_comoving_distance(self) -> np.ndarray[float]:
        """The angular separation scaled by the comoving distance."""
        return np.array(self.ang_sep.value * self.distance)

    @property
    def angular_distance_separation(self) -> np.ndarray[float]:
        """The projected angular size distance of every galaxy to the center."""
        return np.array(
            self.ang_sep.value * (self.distance / (1 + self.group_redshift))
        )

    @property
    def ranks(self) -> np.ndarray[int]:
        """returns the rank of the galaxies with 1 being the closest"""
        return np.argsort(np.argsort(self.angular_separation)) + 1

    def to_data_frame(self) -> pd.DataFrame:
        """Writes a data frame of all the separation metrics."""
        values = {
            f"Sep{self.name}": self.angular_separation,
            f"AngSep{self.name}": self.angular_distance_separation,
            f"CoSep{self.name}": self.projected_comoving_distance,
            f"Rank{self.name}": self.ranks,
        }
        return pd.DataFrame(values)


def add_separation_metrics(
    galaxy_dmu: pd.DataFrame, group_dmu: pd.DataFrame
) -> pd.DataFrame:
    """
    Works out the ranks and projected angular, comoving, and angular size separation of every
    galaxy in the DMU from the BCG, itercen, and flux cen as well as rank
    """
    new_galaxy_dmu = galaxy_dmu.copy()
    groups = np.array(group_dmu["group_id"])
    sep_dfs = []
    for group in groups:
        local_group = group_dmu[group_dmu["group_id"] == group]
        local_group_redshift = np.array(local_group["median_redshift"])
        local_group_galaxies = galaxy_dmu[galaxy_dmu["group_ids"] == group]

        center_bcg = SkyCoord(
            ra=np.array(local_group["bcg_ras"]) * u.deg,
            dec=np.array(local_group["bcg_decs"]) * u.deg,
        )
        center_iter = SkyCoord(
            ra=np.array(local_group["iter_ra"]) * u.deg,
            dec=np.array(local_group["iter_dec"]) * u.deg,
        )
        center_col = SkyCoord(
            ra=np.array(local_group["center_of_light_ras"]) * u.deg,
            dec=np.array(local_group["center_of_light_decs"]) * u.deg,
        )
        galaxy_positions = SkyCoord(
            ra=np.array(local_group_galaxies["RAcen"]) * u.deg,
            dec=np.array(local_group_galaxies["Deccen"]) * u.deg,
        )

        bcg_seps = Separations(
            center_bcg, galaxy_positions, local_group_redshift, "BCG"
        )
        iter_seps = Separations(
            center_iter, galaxy_positions, local_group_redshift, "IterCen"
        )
        col_seps = Separations(
            center_col, galaxy_positions, local_group_redshift, "Cen"
        )

        sep_df = pd.concat(
            [
                bcg_seps.to_data_frame(),
                iter_seps.to_data_frame(),
                col_seps.to_data_frame(),
            ],
            axis=1,
        )
        sep_df["UberID"] = np.array(local_group_galaxies["UberID"])
        sep_dfs.append(sep_df)
    final_separation_df = pd.concat(sep_dfs)
    df_merged = new_galaxy_dmu.merge(final_separation_df, on="UberID", how="left")
    return df_merged.fillna(-999)

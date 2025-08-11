"""
Data structures storing the UCD and Units of all the values.
"""
from dataclasses import dataclass

@dataclass
class Unit:
    """Data structure storing the unit, ucd, and description for the GAMA DMU"""
    unit: str
    ucd: str
    description: str


remames_groups = {
    "group_id": "GroupID",
    "multiplicity": "Nfof",
    "iter_ra": "IterCenRA",
    "iter_dec": "IterCenDec",
    "iter_redshift": "IterCenZ",
    "iter_uber_id": "IterCenUberID",
    "median_redshift": "Zfof",
    "r50": "Rad50",
    "r100": "Rad100",
    "rsimga": "Rad1Sig",
    "velocity_dispersion_gapper": "VelDisp",
    "velocity_dispersion_gaqp_err": "VelDispErr",
    "lum_corrected_mass": "MassProxy",
    "lum_corrected_flux": "TotFluxProxy",
    "bcg_uber_id": "BCGUberID",
    "bcg_ras": "BCGRA",
    "bcg_decs": "BCGDec",
    "bcg_redshifts": "BCGZ",
    "center_of_light_ras": "CenRA",
    "center_of_light_decs": "CenDec",
    "total_absolute_mag": "TotRmag",
    "MassA": "MassA",
    "LumB": "LumB",
    "MassAfunc": "MassAfunc",
    "LumBfunc": "LumBfunc"
}

renames_pairs = {
    'pair_id': "PairID",
    'projected_separation': "SepProjAng",
    'velocity_separation': "SepVel",
    'ra_bar': "RAbar",
    'dec_bar': "Decbar",
    'redshift_bar': "Zbar",
    'total_absolute_mag': "MagSum",
    'uber_id_1': "UberID1",
    'uber_id_2': "UberID2",
}

rename_galaxies = {
    'UberID',
    'RAcen',
    'Deccen',
    'Z',
    'flux_rl',
    'NQ',
    'apparent_mags': Unit("mag", ""),
    'absolute_mags',
    'GroupID': Unit("-", "meta.id", "GroupID for this galaxy in G3CFoFGroup, 0 is un-grouped"),
    'SepBCG': Unit("arcsec", "phys.angSize", "	Projected angular separation of the galaxy to the RA and Dec of the group `BCG` coordinates"),
    'AngSepBCG': Unit("Mpc/h", "phys.size.radius", "	Projected angular size distance (comoving/(1+z)) separation of the galaxy to the RA and Dec of the group `BCG` coordinates"),
    'CoSepBCG': Unit("Mpc/h", "phys.size.radius", "Projected comoving distance separation of the galaxy to the RA and Dec of the group `BCG` coordinates"),
    'RankBCG': Unit("-", "-", "Relative rank of the galaxy to the RA and Dec of the group `BCG` coordinates, where 1 indicates closest and 2 second closest etc"),
    'SepIterCen': Unit("arcsec", "phys.angSize", "Projected angular separation of the galaxy to the RA and Dec of the group `IterCen` coordinates"),
    'AngSepIterCen': Unit("Mpc/h", "phys.size.radius", "Projected angular size distance (comoving/(1+z)) separation of the galaxy to the RA and Dec of the group `IterCen` coordinates"),
    'CoSepIterCen': Unit("Mpc/h", "phys.size.radius", "Projected comoving distance separation of the galaxy to the RA and Dec of the group `IterCen` coordinates"),
    'RankIterCen': Unit("-", "-", "Relative rank of the galaxy to the RA and Dec of the group `IterCen` coordinates, where 1 indicates closest and 2 second closest etc"),
    'SepCen': Unit("arcsec", "phs.angSize", "Projected angular separation of the galaxy to the RA and Dec of the group `Cen` coordinates"),
    'AngSepCen': Unit("Mpc/h", "phys.size.radius", "Projected angular size distance (comoving/(1+z)) separation of the galaxy to the RA and Dec of the group `Cen` coordinates"),
    'CoSepCen': Unit("Mpc/h", "phys.size.radius", "Projected comoving distance separation of the galaxy to the RA and Dec of the group `Cen` coordinates"),
    'RankCen': Unit("-", "-", "Relative rank of the galaxy to the RA and Dec of the group `Cen` coordinates, where 1 indicates closest and 2 second closest etc")
}

group_unit = {
    "GroupID": Unit("-", "-", "Unique group ID"),
    "Nfof": Unit("-", "-", "Group multiplicity"),
    "IterCenRA": Unit("deg", "pos.eq.ra;meta.main", "RA of the iterative central galaxy (J2000)"),
    "IterCenDec": Unit("deg", "pos.eq.dec;meta.main", "Dec of the iterative central galaxy (J2000)"),
    "IterCenZ": Unit("-", "src.redshift", "Redshift of the iterative central galaxy"),
    "IterCenUberID": Unit("-", "meta.id", "Reference ID of the iterative central galaxy (preferred position of the group centre)"),
    "Zfof": Unit("-", "src.redshift;meta.main", "Median redshift of the group"),
    "Rad50": Unit("Mpc/h", "phys.size.radius", "Group radius defined by the 50th percentile group member, based on the projected distance away from IterCenID"),
    "Rad100": Unit("Mpc/h", "phys.size.radius", "Group radius defined by the most distant group member, based on the projected distance away from IterCenID"),
    "Rad1Sig": Unit("Mpc/h", "phys.size.radius", "Group radius defined by the 68th percentile group member, based on the projected distance away from IterCenID"),
    "VelDisp": Unit("km/s", "phys.veloc", "Group velocity dispersion, corrected for the total group velocity dispersion error (VelErr); set to zero if VelErr>VelDispRaw"),
    "VelDispErr": Unit("km/s", "phys.veloc", "Total group velocity dispersion error"),
    "MassProxy": Unit("-", "-", "Proxy for dynamical group mass estimated via M propto Rad50 * VelDisp^2; value listed corresponds to A=1 in Eq. 18 of Robotham et al. (2011), i.e. calibrated mass is given by A*M, with A the scaling factor described in section 4.3 of Robotham et al. (2011)"),
    "TotFluxProxy": Unit("-", "-", "Proxy for group total r-band luminosity down to M_r - 5 log_10 h = -14 in solar luminosities. It is based on the total observed r-band flux and the fraction of light (as estimated from the global galaxy LF) observable at the group redshift; value listed corresponds to B=1 in Eq. 22 of Robotham et al. (2011)"),
    "BCGUberID": Unit("-", "meta.id", "GAMA ID of BCG"),
    "BCGRA": Unit("deg", "pos.eq.ra", "RA of BCG (J2000)"),
    "BCGDec": Unit("deg", "pos.eq.dec", "Dec of BCG (J2000)"),
    "BCGZ": Unit("-", "src.redshift", "Z of BCG"),
    "CenRA": Unit("deg", "pos.eq.ra", "RA of group for a r-band luminosity weighted CoM of the system (J2000)"),
    "CenDec": Unit("deg", "pos.eq.dec", "Dec of group for a r-band luminosity weighted CoM of the system (J2000)"),
    "TotRmag": Unit("mag - 5log(h)", "phys.magAbs", "r-band absolute magnitude of the group, obtained by summing the r-band luminosities of its members"),
    "MassA": Unit("Msun/h", "phys.mass", "Mass proxy times the global A factor required to get a median unbaised halo mass estimate (A=10.0, see Robotham et al. 2011 section 4.3 for details)"),
    "LumB": Unit("Lsun/h^2", "phys.luminosity", "TotFluxProxy times the global B factor required to get a median unbaised r-band luminosity estimate (B=1.04, see Robotham et al. 2011 section 4.4 for details)"),
    "MassAfunc": Unit("Msun/h", "phys.mass", "Mass proxy times the functional A factor which is a function of Nfof and IterCenZ (see Robotham et al. 2011 section 4.3 for details)"),
    "LumBfunc": Unit("Lsun/h^2", "phys.luminosity", "TotFluxProxy times the functional B factor which is a function of Nfof and IterCenZ (see Robotham et al. 2011 section 4.4 for details)"),
}

pair_unit = {
    "PairID": Unit("-", "meta.id;meta.main", "Unique paid ID (Same as GroupID)"),
    "SepProjAng": Unit("kpc/h", "pos.distance", "Physical projected separation between galaxy 1 and galaxy 2"),
    "SepVel": Unit("km/s", "phys.veloc", "Velocity separation of galaxy 1 and galaxy 2"),
    "RAbar": Unit("deg", "pos.eq.ra", "RA of the r-band (petro) flux weighted centre-of-light of the pair"),
    "Decbar": Unit("deg", "pos.eq.dec", "Dec of the r-band (petro) flux weighted centre-of-light of the pair"),
    "Zbar": Unit("-", "src.redshift", "	Redshift of the r-band (petro) flux weighted centre-of-light of the pair"),
    "MagSum": Unit("mag", "phot.mag", "Total r-band (petro) apparent magnitude of the pair"),
    "UberID1": Unit("-", "meta.id", "UberID of galaxy 2 in G3CGal"),
    "UberID2": Unit("-", "meta.id", "UberID of galaxy 2 in G3CGal"),
}
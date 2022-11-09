"""
Provide access to photometric and spectroscopic data for the SDSS Main Galaxy
Sample (MGS) and the quasi-stellar object (QSO) sample observed during
the SDSS-I and SDSS-II Legacy Survey.

This package requires access to two large data files storing data from SDSS
CasJobs queries collecting metadata, photometric data, and spectroscopic
redshift estimates for the MGS and QSO samples, in Feather format.  The path to
a directory holding these files must be specified in a configuration file.  By
default, the configuration file path is:

  ~/.sdss_gal_qso.config

This may be overriden by defining an SDSS_GAL_QSO_CONFIG environment variable
containing the (full) config file path.  E.g., using the bash shell,

  export SDSS_GAL_QSO_CONFIG='/Users/Me/Configs/sdss.config'

If the config file does not exist when this module is loaded, it is created with
default paths for the data files in place, and it will complain if the files are
not found.  Edit the config file to use the correct paths (the user home
directory alias "~/" may be used).

This package also optionally accesses spectral data for catalog objects, 
provided by SDSS in FITS files that are fetched from SDSS.org and copied
locally as needed.  These files are organized in directories according to
SDSS plate number.  The config file must specify the location of a
directory holding the plate directories and FITS files.


2021-01-02:  Created by Tom Loredo, consolidating earlier experiments
2022-11-08:  Major update to use updated SDSS queries for MGS/LRGs/QSOs
"""

import sys
import os
from configparser import ConfigParser

from .catalog import SpecPhotoCatalog


__all__ = ['mgs', 'lrg', 'qso']


# Have not checked earlier versions, but no obvious problems for 3.6.
__minimum_python_version__ = "3.8"


class UnsupportedPythonError(Exception):
    pass


if sys.version_info < tuple((int(val) for val in __minimum_python_version__.split('.'))):
    raise UnsupportedPythonError("sdss_legacy does not support Python < {}".format(__minimum_python_version__))


def find_data(deprecated=False):
    """
    Use the config file to find the data file paths and verify the files exist.

    If the config file doesn't exist, a default config file will be written at
    the default config file location, or the path in the SDSS_GAL_QSO_CONFIG
    environment variable if it is defined.
    """
    default_config_path = '~/.sdss_legacy.config'
    # Paths for tables from deprecated queries, for backwards compatibility:
    default_gal_path = '~/SDSS_Legacy_Data/Spec_Photo_Gal.feather'
    default_qso_old_path = '~/SDSS_Legacy_Data/Spec_Photo_QSO-old.feather'
    # Current table paths:
    default_mgs_path = '~/SDSS_Legacy_Data/Spec_Photo_MGS.feather'
    default_lrg_path = '~/SDSS_Legacy_Data/Spec_Photo_LRG.feather'
    default_qso_path = '~/SDSS_Legacy_Data/Spec_Photo_QSO.feather'
    # path to the folder of "Lite" spectrum FITS files (accum. spectra):
    default_spec_path = '~/SDSS_Legacy_Data/spec_lite/'

    try:
        config_path = os.path.expanduser(os.environ['SDSS_LEGACY_CONFIG'])
    except KeyError:
        config_path = os.path.expanduser(default_config_path)

    # print('Config file path:', config_path)

    if not os.path.exists(config_path):
        config = ConfigParser()
        config['DEFAULT'] = {'Spec_Photo_Gal': default_gal_path,
                             'Spec_Photo_QSO': default_qso_path,
                             'Spec_Lite': default_spec_path}
        with open(config_path, 'w') as cfile:
            config.write(cfile)

    config = ConfigParser()
    config.read(config_path)
    mgs_path = os.path.expanduser(config['DEFAULT']['spec_photo_mgs'])
    lrg_path = os.path.expanduser(config['DEFAULT']['spec_photo_lrg'])
    qso_path = os.path.expanduser(config['DEFAULT']['spec_photo_qso'])
    spec_path = os.path.expanduser(config['DEFAULT']['spec_lite'])
    # Paths for deprecated query data:
    gal_path = os.path.expanduser(config['DEFAULT']['spec_photo_gal'])
    qso_old_path = os.path.expanduser(config['DEFAULT']['spec_photo_qso_old'])

    # print('MGS data path:', mgs_path)
    # print('LRG data path:', lrg_path)
    # print('QSO data path:', qso_path)

    if not os.path.exists(mgs_path) or not os.path.exists(lrg_path) \
            or not os.path.exists(qso_path):
        raise ValueError('MGS, LRG, and/or QSO data file missing!')

    if not os.path.exists(spec_path):
        raise ValueError('Folder for lite spectral FITS files missing!')

    if not deprecated:
        return mgs_path, lrg_path, qso_path, spec_path
    else:
        return gal_path, qso_old_path, spec_path


mgs_path, lrg_path, qso_path, spec_path = find_data()
mgs = SpecPhotoCatalog(mgs_path, spec_path)
lrg = SpecPhotoCatalog(lrg_path, spec_path)
qso = SpecPhotoCatalog(qso_path, spec_path)


def deprecated_catalogs():
    """
    Return catalog instances for galaxy and QSO galaxies from deprecated
    queries.  The deprecated queries included fewer fields/attributes per
    object, and the `gal` sample did not correspond fully to the MGS.
    """
    gal_path, qso_old_path, spec_path = find_data(deprecated=True)

    gal = SpecPhotoCatalog(gal_path, spec_path)
    qso_old = SpecPhotoCatalog(qso_old_path, spec_path)

    return gal, qso_old

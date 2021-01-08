"""
Accessor classes for SDSS Legacy Survey catalog and object data.
"""

from sys import getsizeof
import os
import urllib.request
import webbrowser

import pandas as pd

from astropy import units as u
from astropy.coordinates import SkyCoord

from .sed_sdss import SED_SDSS


class SpecPhotoObject:
    """
    Provide attribute-based access to catalog data for a single object, via
    references to its parent catalog.  Also provide access to the associated
    spectral data.
    """

    # URL for "lite" (accumulated, vs. mult exposure) spectrum FITS files for
    # legacy spectra in DR16.  See bulk data download instructions:
    # https://www.sdss.org/dr16/data_access/bulk/#OpticalSpectraPer-ObjectLiteFiles
    spec_lite_url = 'https://data.sdss.org/sas/dr16/sdss/spectro/redux/26/spectra/lite/'

    attr_map = dict(oclass='class', subclass='subClass',
                    u='modelMag_u', uErr='modelMagErr_u',
                    g='modelMag_g', gErr='modelMagErr_g',
                    r='modelMag_r', rErr='modelMagErr_r',
                    i='modelMag_i', iErr='modelMagErr_i',
                    zmag='modelMag_z', zmagErr='modelMagErr_z',
                    comments='comments_person')

    def __init__(self, id, row, spec_dir):
        """
        Load object data from a Pandas row (series) instance referencing an
        object's catalog data.
        """
        self.id = id  # bestObjID
        self.row = row

        # Spectrum data info:
        self.plate_dir = os.path.join(spec_dir,
                                      '{:0>4d}'.format(row.plate))
        self.sed_fits_name = 'spec-{:0>4d}-{:0>5d}-{:0>4d}.fits'.format(
            row.plate, row.mjd, row.fiberID)

        # Astropy coord instance, using FK5 frame:
        self.coord = SkyCoord(ra=self.ra*u.degree, dec=self.dec*u.degree,
                              frame='fk5')

    def __getattr__(self, name):
        """
        Pass unresolved attribute calls to the row (a series object).  Also,
        map some attributes to new names, to handle collisions (e.g., 'class')
        or for convenience.

        Finally, this method handles lazy creation of the `sed` attribute.
        """
        if name == 'sed':
            sed = self._get_sed()
            setattr(self, 'sed', sed)
            return sed
        elif name in self.attr_map.keys():
            return self.row[self.attr_map[name]]
        else:
            return getattr(self.row, name)

    def _get_sed(self):
        """
        Load spectral energy distribution (SED) data for the object, using an
        existing FITS file if available, otherwise fetching the file from
        SDSS.org.
        """
        if not os.path.isdir(self.plate_dir):
            print('Plate folder does not exist; creating...')
            os.mkdir(self.plate_dir)
        sed_fits_path = os.path.join(self.plate_dir, self.sed_fits_name)
        # print('FITS path:', sed_fits_path)

        if not os.path.exists(sed_fits_path):
            remote = self.spec_lite_url + '{:0>4d}/'.format(self.row.plate) + \
                self.sed_fits_name
            print('Fetching FITS file:', remote)
            with urllib.request.urlopen(remote) as response, \
                    open(sed_fits_path, 'wb') as ffile:
                data = response.read()  # a bytes object
                ffile.write(data)

        return SED_SDSS(sed_fits_path)

    def xp(self, **kwds):
        """
        Show the DR16 Explorer page for the target object in the default web
        browser.

        If keyword arguments are provided, they are passed to the `webbrowser`
        module's `open` function after an Explorer URL argument.
        """
        url = 'https://skyserver.sdss.org/dr16/en/tools/explore/' + \
            'summary.aspx?sid={}&apid='.format(self.specObjID)
        if kwds:
            webbrowser.open(url, **kwds)
        else:
            webbrowser.open_new_tab(url)

    def cutout(self, **kwds):
        """
        Display a cutout of the sky in the vicinity of the target object in the
        default web browser, with nearby specObj targets indicated by boxes.
        """
        # Use default .4" pixel scale, 512x512 px.  For syntax see:
        # http://skyserver.sdss.org/dr12/en/help/docs/api.aspx
        url = 'http://skyserver.sdss.org/dr16/SkyserverWS/ImgCutout/getjpeg?' +\
            'ra={}&dec={}&width=512&height=512&opt=SG'.format(self.ra, self.dec)
        webbrowser.open_new_tab(url)


class SpecPhotoCatalog:
    """
    Provide attribute-based access to the contents of a dataframe holding
    SDSS catalog data from a query of the joined SpecObj and PhotoObj tables,
    and object accessors that reference the catalog data and read spectral
    data from FITS files.
    """

    # Map aliases of data columns to corresponding column names.
    col_map = dict(oclass='class', id='bestObjID', subclass='subClass',
                   u='modelMag_u', uErr='modelMagErr_u',
                   g='modelMag_g', gErr='modelMagErr_g',
                   r='modelMag_r', rErr='modelMagErr_r',
                   i='modelMag_i', iErr='modelMagErr_i',
                   zmag='modelMag_z', zmagErr='modelMagErr_z',
                   comments='comments_person')

    def __init__(self, catalog_path, spec_path):
        """
        Load data from a Feather format tabular datafile at the specified path.
        """
        self.catalog_path = catalog_path
        self.spec_dir = spec_path
        self.df = pd.read_feather(catalog_path)
        self.df.set_index('bestObjID', inplace=True)
        self.cols = tuple(self.df.columns)
        self.shape = self.df.shape
        self.n = self.shape[0]
        self.size = self.df.size  # number of data elements
        self.mem = getsizeof(self.df)  # bytes in memory

    def __getattr__(self, name):
        """
        Pass unresolved attribute calls to the dataframe.  Also, map some
        dataframe columns to new names, to handle collisions (e.g., 'class') or
        for convenience.
        """
        if name in self.col_map.keys():
            return self.df[self.col_map[name]]
        else:
            return getattr(self.df, name)

    def __getitem__(self, id):
        """
        Access object data indexed by bestObjID (which matches the ID in the
        PhotoObj table).
        """
        return SpecPhotoObject(id, self.df.loc[id], self.spec_dir)

    def ith(self, i):
        """
        Access object data for the i'th object in the table.
        """
        id = self.df.index[i]
        return SpecPhotoObject(id, self.df.loc[id], self.spec_dir)

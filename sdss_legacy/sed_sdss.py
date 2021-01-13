""" 
Module providing a container class for SDSS SED data, including support for
Galactic (Milky Way) extinction corrections.

The SED_SDSS class is specifically for SEDs in the DR 12 Main Galaxy Sample
(MGS), comprising targets in the SDSS-I and SDSS-II Legacy Surveys.  Later
surveys may have data in an incompatible format (e.g., the handling of data from
FITS HDU 2 will not work for BOSS survey targets).

The class loads only the coadded spectrum for a source, not per-exposure
spectra.  It needs only the "lite" SED FITS files.  It reads files following
the *spec* data model documented here:

  https://data.sdss.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html


This module was developed using the following conda env:

    $ conda create --name py38astro python=3.8 anaconda
    $ pip install synphot
    $ pip install dust_extinction
    $ pip install dustmaps
    $ pip install astroquery

Dust maps need to be fetched and stored in a suitable location before using
the `dustmaps` module; see:

  Installation — dustmaps v1.0.4 documentation
  https://dustmaps.readthedocs.io/en/latest/installation.html

The SFD and Planck maps are needed.  These were fetched (from the Harvard
Dataverse research repository) by running the following commands in an IPython
session, with `PATH-TO-MAPS` pointing to a convenient location for the map FITS
files:

    from dustmaps.config import config
    config['data_dir'] = 'PATH-TO-MAPS'

    import dustmaps.sfd  # ~134 MB
    dustmaps.sfd.fetch()

    import dustmaps.planck  # ~1.6 GB
    dustmaps.planck.fetch()

These require ~1.7 GB of space.  The Planck map may take a while to download;
the connection ran at 0.3 to 1 MB/s during development.

The Planck map appears to use Planck's 2013 PR1 map:

  Planck PR1 All Sky Maps
  https://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/index.html

See also:

  Foreground maps - Planck Legacy Archive Wiki
  https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/Foreground_maps

Note that `dustmaps` stores the map location in a configuration file in the
user's home directory, so the maps directory path needn't be hard-coded into
scripts and modules.

Created Nov 18, 2020 by Tom Loredo - Refactored from a script
"""
import re
import webbrowser

from numpy import *

from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from dust_extinction.parameter_averages import CCM89, F99
from dustmaps.sfd import SFDQuery
from dustmaps.planck import PlanckQuery


__all__ = ['SED_SDSS']

# Angstrom to microns conversion:
ang2u = 1.e-4

# Dustmaps:
try:
    sfd = SFDQuery()
    planck = PlanckQuery()
except FileNotFoundError:
    raise FileNotFoundError('Requested dustmap(s) are not installed; see README for instructions!')


class SDSSPhotometry:
    """
    Container for SDSS photometry data, including conversion from nanomaggies
    to asinh magnitudes.
    """

    # Constant for magnitude from flux:
    c_mag = -2.5/log(10.)
    # asinh magnitude constants from https://www.sdss.org/dr12/algorithms/magnitudes/:
    b = array([1.4e-10, 0.9e-10, 1.2e-10, 1.8e-10, 7.4e-10])
    m0 = array([24.63, 25.11, 24.80, 24.36, 22.83])
    bb = 2e9*b

    @staticmethod
    def nm2amag(nm):
        """
        Convert flux in nanomaggies, `nm`, to asinh magnitudes.
        """
        SP = SDSSPhotometry
        return SP.c_mag*(arcsinh(nm/SP.bb) + log(SP.b))

    def __init__(self, nm, nm_ivar):
        """
        Store an array of 5 flux values in `nm` and their errors in `nm_ivar`,
        providing both array and attribute access, and conversion to asinh
        magnitudes.
        """
        # TODO:  Compute amag errs from nm errs using delta method
        self.nm = nm
        self.nm_ivar = nm_ivar
        self.amag = SDSSPhotometry.nm2amag(nm)
        self.u = self.amag[0]
        self.g = self.amag[1]
        self.r = self.amag[2]
        self.i = self.amag[3]
        self.z = self.amag[4]


class SpecLine:
    """
    Container for spectral line information from an SDSS SED FITS file.
    """

    def __init__(self, rec):
        """
        Store line data from a FITS data record.
        """
        self.name = rec[3]  # line name
        self.wl = rec[4]  # wavelength
        self.z = rec[5]  # redshift
        self.z_err = rec[6]
        self.sig = rec[7]  # sigma
        self.sig_err = rec[8]
        self.area = rec[9]  # line area
        self.area_err = rec[10]
        self.ew = rec[11]  # EW
        self.ew_err = rec[12]
        self.clev = rec[13]  # continuum level
        self.clev_err = rec[14]

        self.dof = rec[17]  # degrees of freedom
        self.chi2 = rec[18]  # chi**2

    def __str__(self):
        s = 'Line:  {}\n'.format(self.name)
        s += '{:.3f} ang, z = {:.4f} +- {:.3e}\n'.format(
            self.wl, self.z, self.z_err)
        s += 'Area = {:.2e} +- {:.2e}\n'.format(self.area, self.area_err)
        s += 'Chi**2 = {:.2f} for {:.1f} dof'.format(self.chi2, self.dof)
        return s


class SED_SDSS:
    """
    Container class for spectral energy distribution (SED) data from the
    Sloan Digital Sky Survey (SDSS).

    This is for data from FITS files with the spec data model, i.e., the
    coadded spectrum and (optionally) frames for individual exposures (omitted
    in the "lite" version).  See:
    https://data.sdss.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html
    """

    # Regexp to get FITS file basename.  Spec file names are of the form
    # spec-PLATE-MJD-FIBER.fits.
    ffbase_re = re.compile('.*spec-(?P<base>.*).fits')

    def __init__(self, fits_path):
        """
        Load SED data from an SDSS FITS file.
        """
        self.fits_path = fits_path
        self.fits_base = self.ffbase_re.match(fits_path).group('base')
        with fits.open(fits_path) as hdus:
            hdu0 = hdus[0]  # FITS PrimaryHDU
            self.hdu0 = hdu0
            # The SpecID is a hash of identifying data for the object and the
            # pipeline software producing the FITS data; see:
            # https://www.sdss.org/dr12/spectro/spectro_basics/
            self.spec_id = hdu0.header['spec_id'].strip()  # unique object ID
            self.pid = hdu0.header['plateid']  # plate ID
            self.fid = hdu0.header['fiberid']  # fiber number, 1 - 640
            self.MJD = hdu0.header['MJD']

            # *** What are these coords?
            # self.RA = hdu0.header['RADEG']
            # self.dec = hdu0.header['DECDEG']

            # Object coords:
            self.equinox = hdu0.header['equinox']
            if self.equinox != 2000.:
                raise ValueError('Equinox should be 2000!')
            self.frame = hdu0.header['radecsys']
            if self.frame != 'FK5':
                raise ValueError('Coordinate frame should be FK5!')
            self.ra = hdu0.header['plug_RA']  # RA in J2000 deg
            self.dec = hdu0.header['plug_dec']  # dec in J2000 deg
            # Astropy coord using FK5 frame:
            self.coord = SkyCoord(ra=self.ra*u.degree, dec=self.dec*u.degree,
                                  frame='fk5')

            # For info about the relationship b/t the ICRS system and the J2000
            # frame, see: https://docs.astropy.org/en/stable/coordinates/definitions.html

            self.funit = hdu0.header['bunit']
            self.qual = hdu0.header['quality']  # night/image quality
            self.nexp = hdu0.header['nexp']  # number of exposures

            if hdu0.header['WFitType'] != 'LOG-LINEAR':
                raise ValueError('{}: FITS file wavelength scale not '
                                 'log-linear!'.format(self.fits_base))
            if hdu0.header['NWOrder'] != 2:
                raise ValueError('{}: FITS file wavelength scale order should '
                                 'be 2!'.format(self.fits_base))
            # Wavelength scale coefficients:
            self.wc0 = hdu0.header['Coeff0']  # Central wavelength (log10) of first pixel
            self.wc1 = hdu0.header['Coeff1']  # log10 dispersion per pixel

            # HDU 1 has the coadded spectrum.
            self.loglam = hdus[1].data['loglam']  # lab-frame log wavelength/ang
            self.lam = 10.**self.loglam  # lab-frame wavelength, angstroms
            self.flux = hdus[1].data['flux']      # SED in flux units
            self.ivar = hdus[1].data['ivar']      # inverse variance
            self.and_mask = hdus[1].data['and_mask']
            self.or_mask = hdus[1].data['or_mask']
            self.wdisp = hdus[1].data['wdisp']      # wavelength dispertion (LSF?)
            self.sky = hdus[1].data['sky']      # inverse variance
            # Model fit used for classification and redshift:
            self.model = hdus[1].data['model']

            # HDU 2 has specObj data for the plate and object.
            # (This would instead need to read spAll data for BOSS spectra.)
            self.specObjReader(hdus[2])

            # # Redshift estimates; for non-QSO targets, use z_nq; see:
            # # https://www.sdss.org/dr13/spectro/caveats/
            # self.z = h2d['z'][0]  # "best" redshift
            # self.z_err = h2d['z_err'][0]  # z error
            # self.z_nq = h2d['z_NoQSO'][0]  # z est excluding QSO templates
            # self.z_nq_err = h2d['z_err_NoQSO'][0]
            # self.z_warn = h2d['zWarning_NoQSO'][0]  # z est excluding QSO templates

            # # Target bitmask names, omitting "Target".
            # self.legacy_1 = h2d['legacy_target1']
            # self.legacy_2 = h2d['legacy_target2']
            # # BOSS targets are *mostly* LRGs and quasars.
            # self.boss_1 = h2d['boss_target1']
            # self.boss_anc1 = h2d['ancillary_target1']
            # self.boss_anc2 = h2d['ancillary_target2']

            # Main sample galaxy indicator; see:
            # https://www.sdss.org/dr12/spectro/targets/
            self.main_samp = (self.legacy_1 & (64 | 128 | 256)) > 0

            # HDU 3 has spZline data from Spectro-1D emission line fits.
            # For line IDs, see:
            # https://data.sdss.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/PLATE4/RUN1D/spZline.html
            self.lines = [SpecLine(rec) for rec in hdus[3].data]

    def hdus(self):
        """
        Opens the FITS file for the SED and returns a FITS HDU list.
        """
        return fits.open(self.fits_path)

    def specObjReader(self, hdu):
        """
        Read data from a FITS HDU using the specObj data model.  If `obj`
        is provided, some of the read values will be inserted into the
        attribute namespace of `obj`.

        This gathers data from HDU 2 in the FITS file, which records a row of
        data using the *specObj* data model for objects in the MGS (note that
        BOSS uses a different format for HDU 2 data).  The specObj data model is
        described here:

        https://data.sdss.org/datamodel/files/SPECTRO_REDUX/specObj.html
        """
        data = hdu.data

        # Spectroscopic ID:
        specobjid = data['specObjID'][0].strip()
        if self.spec_id != specobjid:
            print('{}: Mismatch bewteen spec_id and specObjID'
                  '....'.format(self.fits_base))

        # SDSS photometric ObjID is an int array of 5 values:
        # RUN, RERUN, CAMCOL, FIELD, ID)
        self.objid = data['ObjID'][0]
        self.run, self.rerun, self.camcol, self.field, self.id = self.objid

        # Related photometric IDs:
        self.fluxobjid = data['fluxObjID'][0]
        # "best" is the recommended ID for joining spec and photom data.
        self.bestobjid = data['bestObjID'][0]
        self.targetobjid = data['targetObjID'][0]

        self.survey = data['survey'][0]
        self.program = data['programName'][0]
        self.pquality = data['plateQuality'][0]

        ra = data['plug_RA'][0]  # RA in J2000 deg
        dec = data['plug_dec'][0]  # dec in J2000 deg
        # Double-check that we got the right coords:
        if ra != self.ra or dec != self.dec:
            print('{}: Mismatch bewteen spec and specObjID RA, dec'
                  '....'.format(self.fits_base))

        # Spectroscopic class'n:
        self.sclass = data['class'][0]  # ".class" is reserved
        self.subclass = data['subclass'][0]
        self.is_star = self.sclass == 'STAR'
        self.is_gal = self.sclass == 'GALAXY'
        self.is_qso = self.sclass == 'QSO'

        # Image class'n (type, or objType in DR8):
        self.itype = data['sourceType'][0]
        # itype indicates why the object was targeted; sclass is a more reliable
        # classification; see the DR12 FAQ:
        # https://www.sdss.org/dr12/help/faq/#spectype2

        # Observed spectral range (ang):
        self.wmin = data['wavemin'][0]
        self.wmax = data['wavemax'][0]

        # Redshift estimates; for BOSS non-QSO targets, use z_nq; see:
        # https://www.sdss.org/dr13/spectro/caveats/#noqso
        self.z = data['z'][0]  # "best" redshift
        self.z_err = data['z_err'][0]

        # # NoQSO results are for survey = 'boss' spectra only.
        # self.z_nq = data['z_NoQSO'][0]  # z est excluding QSO templates
        # self.z_nq_err = data['z_err_NoQSO'][0]
        # # Warning bitmask; 0 = OK:
        # self.z_warn = data['zWarning_NoQSO'][0]  # z est excluding QSO templates

        self.vdisp = data['vdisp'][0]  # velocity dispersion (km/s)
        self.vdisp_err = data['vdisp_err'][0]

        # Target bitmask names, omitting "Target".
        self.legacy_1 = data['legacy_target1'][0]
        self.legacy_2 = data['legacy_target2'][0]
        # BOSS targets are *mostly* LRGs and quasars; expect none in MGS.
        self.boss_1 = data['boss_target1'][0]
        self.boss_anc1 = data['ancillary_target1'][0]
        self.boss_anc2 = data['ancillary_target2'][0]

        # Photometric info:
        # Spectrum projected onto ugriz filters (nanomaggies); 5-vector:
        self.sphotom = SDSSPhotometry(data['spectroflux'][0],
                                      data['spectroflux_ivar'][0])

        # Best-fit template photometry:
        self.tphotom = SDSSPhotometry(data['spectrosynflux'][0],
                                      data['spectrosynflux_ivar'][0])

        # Calibration photometry:
        self.cphotom = SDSSPhotometry(data['calibflux'][0],
                                      data['calibflux_ivar'][0])

        # Human visual inspection results:
        # Redshift:
        self.hz = data['z_person'][0]
        # Redshift conficence (0=not inspected or no confidence, 1,2=low
        # confidence, 3,4=high confidence):
        self.hz = data['z_person'][0]
        # Class (0=not inspected or unknown, 1=star, 2=narrow emission-line
        # galaxy, 3=QSO, 4=galaxy):
        self.hclass = data['class_person'][0]
        # Comments:
        self.comments = data['comments_person'][0]

    def dereddened(self, Rv=3.1, curves='CCM89', map='SFD', factor=False):
        """
        Adjust the flux for Galactic extinction, using the provided reddening
        parameter and choice of extinction curve family ('CCM89' or 'F99') and
        dust map ('SFD' or 'Planck').
        """
        if map.lower() == 'sfd':
            Ebv = sfd(self.coord)
        elif map.lower() == 'planck':
            Ebv = planck(self.coord)
        else:
            raise ValueError('Invalid dust map selection!')

        if curves.lower() == 'ccm89':
            ext = CCM89(Rv=Rv)
        elif curves.lower() == 'f99':
            ext = F99(Rv=Rv)
        else:
            raise ValueError('Invalid extinction curve selection!')

        # Note that ext. curves use wavelength in microns
        efactor = ext.extinguish(ang2u*self.lam*u.micron, Ebv=Ebv)
        flux = self.flux/efactor
        ivar = self.ivar*efactor**2

        if factor:
            return flux, ivar, efactor
        else:
            return flux, ivar

    def xp(self, **kwds):
        """
        Show the DR16 Explorer page for the target object in the default web
        browser.

        If keyword arguments are provided, they are passed to the `webbrowser`
        module's `open` function after an Explorer URL argument.
        """
        url = 'https://skyserver.sdss.org/dr16/en/tools/explore/' + \
            'summary.aspx?sid={}&apid='.format(self.spec_id)
        if kwds:
            webbrowser.open(url, **kwds)
        else:
            webbrowser.open_new_tab(url)

    def cutout(self, **kwds):
        """
        Display a cutout of the sky in the vicinity of the target object in the
        default web browser, with specObj targets indicated by boxes.
        """
        # Use default .4" pixel scale, 512x512 px.  For syntax see:
        # http://skyserver.sdss.org/dr12/en/help/docs/api.aspx
        url = 'http://skyserver.sdss.org/dr16/SkyserverWS/ImgCutout/getjpeg?' +\
            'ra={}&dec={}&width=512&height=512&opt=SG'.format(self.ra, self.dec)
        webbrowser.open_new_tab(url)

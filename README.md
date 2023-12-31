# `sdss_legacy`: Access photometric and spectroscopic data from the SDSS-I/SDSS-II Legacy Survey

The `sdss_legacy` package provides access to photometric and spectroscopic data (including SEDs) for the SDSS Main Galaxy Sample (MGS), the luminous red galaxy (LRG) sample, and the quasi-stellar object (QSO) sample observed during the SDSS-I/SDSS-II Legacy Survey.  The MGS sample comprises 678,239 objects (mainly regular galaxies, but including some non-QSO AGNs, e.g., Seyferts, and a subset of LRGs). The LRG sample comprises 187,423 objects (some in the MGS sample). The QSO sample comprises 74,389 objects.

Specifically, this package accesses the Legacy data as provided by SDSS DR17 (the last data release of SDSS IV). For basic information about the samples, see: [Legacy Survey Target Selection | SDSS DR17](https://www.sdss.org/dr17/algorithms/legacy_target_selection/). Note that although these samples are "legacy" samples defined by observations taken during SDSS I/II, pipeline revisions lead to small changes in object properties and thus sample selections across data releases.

The data tables used by `sdss_legacy` were created using CasJobs SQL queries building an inner join between the SDSS `SpecObj` and `PhotoObj` tables, and a left outer join with the `Photoz` table (providing SDSS $k$NN photo-z estimates for a subset of the objects). The resulting samples are about 1% smaller than samples with the same `WHERE` selection clause operating on just the `SpecObj` table. We have not explored how omitting the lost samples (lacking `PhotoObj` entries) may introduce or reduce bias in the samples for certain purposes. The queries are specified to include morphological summaries from the photometric pipeline, as well as various magnitudes (e.g., model-based, PSF, Petrosian, and fiber), and basic spectroscopic summaries (e.g., redshift, velocity dispersion). The queries used to build the samples appear at the end of this document.



---

---

## Installation and setup

### Usage/development environment

This package is being developed using the following `conda` environment; use the package in a similar environment:
```bash
$ conda create --name py310astro -c conda-forge -c defaults ipython scipy matplotlib \
  h5py beautifulsoup4 html5lib bleach pandas sortedcontainers \
  pytz setuptools mpmath bottleneck jplephem asdf pyarrow colorcet
$ conda activate py310astro
$ pip install synphot
$ pip install dust_extinction
$ pip install dustmaps
$ pip install astroquery
$ pip install speclite
```

(The `conda create` command defines a basic environment for AstroPy, used for FITS file access; this environment definition supplements the basic AstroPy requirements with `ipython`.)

The `sdss_legacy` package may be installed using `pip` and the path to the distribution (i.e., the directory containing this README and the `setup.py` file). E.g., from within a directory containing the distribution:

```bash
pip install sdss_legacy
```

If you plan to modify the package, consider installing in developer mode:

```bash
pip install -e sdss_legacy
```

This will put links to the modules in the distribution in the current Python environments `site-packages` directory, rather than copying the modules there. With this setup, changes to the modules take effect immediately on the next import; there is no need to re-install the modified package.

### Dust maps

For Galactic (Milky Way) extinction corrections, `sdss_legacy` uses the [dust_extinction](https://dust-extinction.readthedocs.io/en/stable/) and [dustmaps](https://dustmaps.readthedocs.io/en/latest/index.html) packages. Dust maps need to be fetched and stored in a suitable location before using the `dustmaps` package (and thus before using `sdss_legacy`). For detailed instructions, see:

> [Installation — dustmaps documentation](https://dustmaps.readthedocs.io/en/latest/installation.html)

Specifically for `sdss_legacy`, the SFD and Planck maps are needed.  These may be fetched (from the Harvard Dataverse research repository) by running the following commands in an IPython session, with `PATH-TO-MAPS` pointing to a convenient location for the map FITS
files:
```python
from dustmaps.config import config
config['data_dir'] = 'PATH-TO-MAPS'

import dustmaps.sfd  # ~134 MB
dustmaps.sfd.fetch()

import dustmaps.planck  # ~1.6 GB
dustmaps.planck.fetch()
```
These require ~1.7 GB of space.  The Planck map may take a while to download;
the connection ran at 0.3 to 1 MB/s during development.

### Object catalog data storage and configuration file

This package requires access to two large data files storing data from SDSS CasJobs SQL queries collecting metadata, photometric data, and spectroscopic redshift estimates for the MGS and QSO samples, in Feather format.  The path to a directory holding these files must be specified in a configuration file.  By default, the configuration file path is:
```
~/.sdss_legacy.config
```

This may be overriden by defining an `SDSS_LEGACY_CONFIG` environment variable containing the (full) config file path.  E.g., using the bash shell,
```bash
export SDSS_LEGACY_CONFIG='/Users/Me/Configs/sdss.config'
```
If the config file does not exist when the package is loaded, it is created with default paths for the data files defined, and Python will complain if the files are not found there when requested.  Edit the config file to use the correct paths (the user home directory alias "~/" may be used in path names), or put the needed Feather-format files (and, optionally, FITS spectral data files) here:

* `~/SDSS_Legacy_Data/` for Feather files with SDSS catalog data
* `~/SDSS_Legacy_Data/spec_lite` for spectral data in FITS files

A good way to set up a config file with the proper format is to install the package, and then just import it in a `python` or `ipython` session:

```python
import sdss_legacy
```

This will raise an exception about missing data files (unless you happened to install them at the default locations), but in the meantime it will create a default config file in the default location. Adjust it as necessary.

For Cornell PZ team users, the data tables are pre-installed on the Tabulator machine:

* Feather files in `/data/photoz/sql_queries`
* FITS spectral data files in `/data/photoz/spec_lite`

Store the files in a convenient location, and modify the config file to point to them. The Feather files were derived from FITS files produced by SQL queries submitted to CasJobs; they are smaller than the FITS files, and they can be loaded much more quickly (in about a second on my i7 Mac Mini). The package loads them as Pandas dataframes; the table dataframe sizes are (rows, memory):

* MGS: 678239 rows,
  624.8 MB
* LRG: 187423 rows, 171.2 MB
* QSO: 74389 rows,
  69.3 MB

Note that there is overlap between the tables: entries with the same `specObjID` identifier appear in multiple tables.

In the future we will likely move the tables to a research data repository such as Zenodo or Dataverse.


### Spectral data FITS file storage

This package also optionally accesses spectral data for catalog objects, provided by SDSS in FITS files that may be fetched from SDSS.org and copied locally as needed.  These files are organized in directories according to SDSS plate number.  The config file must specify the location of a directory holding the plate directories and FITS files. This directory can be empty; it will get populated by FITS files as they are requested. Alternatively, you may download them in bulk.

Information about bulk data downloads is available here: [Bulk Data Downloads | SDSS DR17](https://www.sdss.org/dr17/data_access/bulk/). The package accesses the data via HTTP requests to this SDSS archive:

> https://data.sdss.org/sas/dr17/sdss/spectro/redux/26/spectra/lite/

Note that this is the "Lite" archive, providing only the accumulated spectrum for an object, not the spectra in the multiple exposures used to make the coadded spectrum.


---

---



## Usage

Beware, in this documentation, the potential for confusion between "object" in the sense of a Python object (or instance), and "object" in the sense of an astronomical object (e.g., a galaxy or QSO).

Import the `mgs`, `lrg`,  and `qso` `SpecPhotoCatalog` instances ("catalogs") from the package:

```python
from sdss_legacy import mgs, lrg, qso
```

The import will load the data tables. The `mgs`, `lrg`,  and `qso` catalogs are wrappers around Pandas dataframes. You may directly access the dataframes via `mgs.df`, `lrg.df`, and `qso.df`. The easiest way to access data columns is as attributes of the catalog objects. E.g., 

```python
mgs.z
```

returns all of the redshifts in the MGS sample (as a Pandas `Series` object).

See the SQL queries below for the available attributes.

Some attributes have been *mapped to alternative names*, to avoid namespace collisions (e.g., with Python's `class` attribute for class instances), or for convenience. For example, the SDSS `class` attribute can be accessed by using `oclass` ("object class") as its identifier:

```python
qso.oclass
```

The Python token `id` maps to the SDSS  `bestObjID` (the index attribute; see below), and `comments` maps to the SDSS  `comments_person`. There are also simplified identifiers for model magnitudes (the recommended magnitude type for galaxies), e.g., `mgs.u` for the *u* magnitudes (and similarly for *g*, *r*, and *i*), and `mgs.zmag` for the *z* magnitudes (distinguishing it from the redshift attribute, `z`); these may change. The associated uncertainties are `mgs.uErr` (etc.) and `mgs.zmagErr` (note that this violates NumPy naming conventions, instead aiming to resemble the SDSS attribute names).

Of course, the underlying dataframe can be used to access columns by SDSS name using strings, e.g., `mgs.df['class']` for the column of galaxy classes.

The catalogs are indexed by the `SpecObj`table's  `bestObjID` identifier (henceforth "ID"), which is the same as the `PhotoObj` `objID` (and different from the `SpecObj` `specObjID`). The name of this attribute reflects that it accesses photometric data for the object in the photometric catalog that "best" matches the location of the fiber used to obtain a spectrum. 

Accessing `mgs.index` or `qso.index` returns the dataframe index (as a Pandas index object).



### Accessing catalog data for an object

Accessing data for an object via its ID produces a `SpecPhotoObject` instance providing access to the object's data. E.g., for the first galaxy in the MGS table, an object accessing the galaxy data may be obtained via:

```python
obj = mgs[1237667323264893037]
```

For convenience, you can access the data for the $i$'th object in a catalog via the `ith()` method (note that this is a method, used with function call syntax, not via `ith[]`, and thus currently without an iterator interface). So data for the same galaxy may be accessed via:

```python
obj2 = mgs.ith(0)
```

These instances reference data in the tables; they don't copy the data (to help constrain memory use).

The object instances are wrappers around Pandas `Series` objects containing labeled data for a row in a table. The data may be accessed using the same attribute names available for columns in the table (as listed in the SQL queries below, and adopting the aliases mentioned above).

An object instance also has a `coord` attribute holding an instance of the Astropy `Coord` object for sky coordinates. This is derived from the J2000 FK5 RA and declination held in the `ra` and `dec` attributes (in degrees). The `coord` attribute is used internally to estimate Galactic extinction in the direction of the object.



### Quick looks at object data

Two convenience functions let you display summaries of an object's data in your default web browser.

```python
obj.xp()
```

will open a new tab displaying the [SDSS DR17 Object Explorer](http://skyserver.sdss.org/dr17/en/tools/explore/Summary.aspx?) page for the object, showing an image thumbnail, a plot of the spectrum, and some key summary data (mnemonic: "xp" = "ExPlore").

```python
obj.cutout()
```

will open a new tab displaying a cutout of the sky in the vicinity of the object (at higher resolution than the Explorer thumbnail), with nearby `specObj` targets indicated by boxes.



### SED data access

A `SpecPhotoObject` instance has a `sed` attribute providing access to spectral data measuring the spectral energy distribution (SED) via an instance of the `SED_SDSS` class. This attribute is defined "lazily;" it is not defined until it is used, so the spectrum FITS file is neither fetched from SDSS.org (if not already downloaded) nor loaded if the user does not access spectral data.

The `sed` instance provides attribute access to a large number of quantities loaded from the spectrum FITS file. Some of the quantities are redundant with table data, but possibly useful for cross-checking (the `SED_SDSS` class is designed to be usable separately from the catalogs, in which case many attributes may be more useful). The spectrum itself, with measurement errors, is available via:

```python
sed.wmin, sed.wmax  # wavelength range (anstroms)
sed.loglam  # array of observer-frame log(wavelength/1 angstrom)
sed.lam  # array of observer-frame wavelength (angstrom)
sed.flux  # array with observer-frame SED in flux units (1E-17 erg/cm^2/s/Ang)
sed.ivar  # array of inverse variances for flux measurement errors
```

Note that the SDSS processing makes some elements of `sed.flux` negative for some objects.

The SDSS pipeline performs a number of analyses of the spectrum, with results available in the FITS file, many accessible via attributes, including the following:

* Redshift (`z` and `z_err`) and velocity dispersion (`vdisp` and `vdisp_err`) are available (duplicating catalog data).
* `lines` is a list of `SpecLine` instances containing results of fits of spectral lines, including these attributes:
  * `name` — line name
  * `wl` — wavelength (angstroms)
  * `z`, `z_err` — line redshift and uncertainty
  * `sig` — line width as standard deviation
  * `area` — area in the line
  * `ew`, `ew_err` — equivalent width and uncertainty (angstroms)
* 3 types of synthetic photometry results (in nanomaggies): `sphotom` ("spectroflux") projecting the SED onto $ugriz$ filters; `tphotom` ("spectrosynflux") with photometry from a best-fit template, and `cphotom` ("calibflux") with "calibration photometry." These are `SDSSPhotometry` instances with attributes providing fluxes (`nm`) and inverse variances (`nm_ivar`) as 5-element arrays (with *ugriz* data). There are also `u`, `g`, `r`, `i`, and `z` attributes with magnitudes corresponding to the flux measurements.

Quite a few other quantities are loaded from the FITS file; see the `sed_sdss.py` source file for details.



### Galactic extinction corrections for SEDs

The `SED_SDSS` class provides capability for correcting a spectrum for Galactic extinction. Dust and gas in the Galaxy attenuates light from extragalactic objects in a direction- and wavelength-dependent manner. The behavior has been measured, and the attenuation for a particular object can largely be described with two empirically identified parameters:

* **Monochromatic attenuation:** Extinction *dims* a spectrum. $A_V$ is the *attenuation* due to extinction at the effective wavelength of the $V$ band, in magnitudes (positive changes in magnitude correspond to dimming).

* **Color excess:** Extinction *reddens* a spectrum: the attenuation tends to be stronger at short (blue) wavelengths than at long (red) wavelengths, making objects appear redder than they actually are. The *color excess* measures how extinction distorts the *shape* of a spectrum. The conventional measure is the (positive) color excess induced in the $B-V$ color; this is denoted $E(B-V)$ or $E_{B-V}$. Extinction changes the $B-V$ color to $(B+A_B) - (V+A_V)$, so  $E_{B-V} = A_B - A_V$.

These two quantities are tied together: along a line of sight through dust of homogeneous composition, the more dust one looks through (i.e., the greater the distance), the greater the attenuation, and the greater the reddening. The relationship between the two depends on the composition of dust along the line of sight, which determines the **extinction curve**: the attenuation as a function of wavelength. Empirically, curves of extinction vs. wavelength for different dust compositions do not cross and are nearly monotonic in wavelength. They thus can be labeled by a single parameter, conventionally taken to be the *reddening parameter*, $R_V$, with $R_V \equiv A_V/E_{B-V}$. A typical value is $R_V = 3.1$, i.e., dust producing a color excess of about a third of of the $V$ band attenuation. (Some studies find that extinction laws are better described using two parameters to specify the shape.)

Astronomers so far have not produced definitive measurements and parameterizations of Galactic extinction. `sdss_legacy` users have to make two choices in order to correct a spectrum for extinction, specifying how to handle direction- and wavelength-dependence of extinction.

* **Direction dependence:** A *dust map* provides an estimate of the amount of dust along a line of sight out of the Galaxy in a particular direction, measured by the color excess, $E_{B-V}$. `SED_SDSS` instances support two choices of dust maps: the Schlegel, Finkbeiner & Davis (1998, 2011) or **SFD** map, and the map measured by the **Planck** mission.
* **Wavelength dependence:**  `SED_SDSS` instances support two choices of *extinction curve families*, each using $R_V$ to parameterize the curves within the family: the Cardelli, Clayton, and Mathis (1989) **CCM89** family, and the Fitzpatrick (1999, 2004) **F99** family. 

Additionally, one must specify a value for the $R_V$ parameter identifying which curve to use in a particular family.

`SED_SDSS` instances provide a `dereddened(Rv=3.1, curves='CCM89', map='SFD')` method that returns a copy of the `flux` and `ivar` data, corrected for Galactic extinction using the specified choices of $R_V$, curve family, and map. The adjusted flux is the measured flux divided by a wavelength-dependent *extinguishing factor* (or *extinction factor*). The adjusted inverse variance is the raw value multiplied by the extinguishing factor squared (this may not be accurate, depending on how much of the measurement error is due to operations like background subtraction). Adding a `factor=True` argument will produce a third return value, the array of extinguishing factors (for each wavelength).



## Example script

Copy `ExtinctionExpt.py` from the `scripts` folder to a location outside of the distribution, and execute it with `ipython -i ExtinctionExpt.py` (the `-i` option makes iPython stay interactive after executing the script, leaving plot windows open). This scripts loads the catalogs, collects data for 10 MGS galaxies, 10 LRG galaxies, and 10 QSOs, and then makes some plots:

* A plot showing spectra of 10 MGS galaxies on a common set of axes.
* 3 plots (1 for an MGS galaxy, 1 for an LRG galaxy, 1 for a QSO), each with 4 panels showing extinction-corrected SEDs for the possible choices of `(curves, map)`. In each panel, a gray line shows the adjusted spectrum, and vertical error bars dip from the curve to the unadjusted spectrum (zoom in to see the error bars more clearly). A dark red curve (against the right axis) shows the extinction factor (the adjusted SED is the measured SED *divided* by this factor).

The plots suggest that the choice of map has more impact than the choice of extinction curve family, but this is only from a cursory exploration.

At the interactive prompt (after the plots appear), the variables `m`, `l`, and `q` reference the catalog objects for the galaxies used in the plots (from the MGS, LRG, and QSO tables, resp.). You can use these to experiment with attribute access and the `xp()` and `cutout()` methods.



---

---

## CasJobs SQL queries

The data tables used by `sdss_legacy` were generated by the following set of three queries (they have identical `SELECT` content, but different `WHERE` clauses and output table names). Use these queries to identify object attribute names (other than the renamed attributes identified above).

```sql
SELECT
    -- specObj metadata:
    CONVERT(bigint, so.specObjID) as specObjID,
    so.fluxObjId, so.bestObjID,
    so.sciencePrimary, so.sdssPrimary, so.legacyPrimary,
    so.firstRelease, so.survey,
    so.programname, so.mjd,
    so.plate, so.fiberID, so.run2d,
    so.ra, so.dec, so.cx, so.cy, so.cz,
    -- specObj results:
    so.z, so.zErr,
    so.class, so.subClass, so.rChi2, so.DOF, so.rChi2Diff,
    so.z_person, so.class_person, so.comments_person,
    so.velDisp, so.velDispErr,
    so.waveMin, so.waveMax, so.wCoverage,
    -- Subsample membership:
    CONVERT(bit, (so.legacy_target1 & 32)/32) as is_lrg,  -- LRG (all criteria)
    CONVERT(bit, (so.legacy_target1 & 64)/64) as is_gal,  -- regular or LSB galaxy (mostly)
    CONVERT(bit, (so.legacy_target1 & 128)/128) as is_lsb, -- LSB galaxy
    CONVERT(bit, (so.legacy_target1 & 256)/256) as is_bc,  -- bright-core
    -- PSF mag:
    po.psfMag_r, po.psfMagErr_r, 
    -- modelMag uses best of deV or exp profiles (max like in r band):
    po.modelMag_u, modelMagErr_u,
    po.modelMag_g, modelMagErr_g, 
    po.modelMag_r, modelMagErr_r,
    po.modelMag_i, modelMagErr_i, 
    po.modelMag_z, modelMagErr_z,
    -- deV profile params:
    po.deVMag_u, po.deVMagErr_u, po.deVRad_u, po.deVAB_u,
    po.deVMag_g, po.deVMagErr_g, po.deVRad_g, po.deVAB_g,
    po.deVMag_r, po.deVMagErr_r, po.deVRad_r, po.deVAB_r,
    po.deVMag_i, po.deVMagErr_i, po.deVRad_i, po.deVAB_i,
    po.deVMag_z, po.deVMagErr_z, po.deVRad_z, po.deVAB_z,
    -- exp profile params:
    po.expMag_u, po.expMagErr_u, po.expRad_u, po.expAB_u,
    po.expMag_g, po.expMagErr_g, po.expRad_g, po.expAB_g,
    po.expMag_r, po.expMagErr_r, po.expRad_r, po.expAB_r,
    po.expMag_i, po.expMagErr_i, po.expRad_i, po.expAB_i,
    po.expMag_z, po.expMagErr_z, po.expRad_z, po.expAB_z,
    -- r-band model likelihoods:
    po.lnLDeV_r, po.lnLExp_r,
    -- petroMag can be used to estimate surface brightness:
    po.petroMag_u, po.petroMagErr_u, po.petroRad_u, po.petroR50_u, po.petroR90_u,
    po.petroMag_g, po.petroMagErr_g, po.petroRad_g, po.petroR50_g, po.petroR90_g,
    po.petroMag_r, po.petroMagErr_r, po.petroRad_r, po.petroR50_r, po.petroR90_r,
    po.petroMag_i, po.petroMagErr_i, po.petroRad_i, po.petroR50_i, po.petroR90_i,
    po.petroMag_z, po.petroMagErr_z, po.petroRad_z, po.petroR50_z, po.petroR90_z,
    -- fiber mags, 3" for SDSS spectrograph, 2" for BOSS:
    po.fiberMag_u, fiberMagErr_u, po.fiberMag_g, fiberMagErr_g, 
    po.fiberMag_r, fiberMagErr_r, po.fiberMag_i, fiberMagErr_i, 
    po.fiberMag_z, fiberMagErr_z,
    po.fiber2Mag_u, fiber2MagErr_u, po.fiber2Mag_g, fiber2MagErr_g, 
    po.fiber2Mag_r, fiber2MagErr_r, po.fiber2Mag_i, fiber2MagErr_i, 
    po.fiber2Mag_z, fiber2MagErr_z,
    -- photo extinction corrections:
    po.extinction_u, po.extinction_g, po.extinction_r,
    po.extinction_i, po.extinction_z,
    -- spec synthetic mags (direct, not template-based; nanomaggies):
    so.spectroFlux_u, so.spectroFluxIvar_u, 
    so.spectroFlux_g, so.spectroFluxIvar_g, 
    so.spectroFlux_r, so.spectroFluxIvar_r, 
    so.spectroFlux_i, so.spectroFluxIvar_i, 
    so.spectroFlux_z, so.spectroFluxIvar_z,
    -- photo-z:
    ISNULL(pz.z, -1) as z_p,
    ISNULL(pz.zErr, -1) as z_p_err,
    ISNULL(pz.nnCount, -1) as z_p_nn,
    ISNULL(pz.photoErrorClass, -1) as z_p_err_class

-- For MGS:
INTO mydb.Spec_Photo_MGS
FROM SpecObj AS so
    JOIN PhotoObj AS po ON po.objID = so.bestObjID
    LEFT OUTER JOIN Photoz pz ON (pz.objID = so.bestObjID)
WHERE
    so.legacyPrimary = 1 AND so.zWarning = 0 AND so.class = 'GALAXY' -- gals w/ good spec
    AND (so.legacy_target1 & (64 | 128 | 256)) > 0  -- MGS

-- For LRGs:
INTO mydb.Spec_Photo_LRG
FROM SpecObj AS so
    JOIN PhotoObj AS po ON po.objID = so.bestObjID
    LEFT OUTER JOIN Photoz pz ON (pz.objID = so.bestObjID)
WHERE
    so.legacyPrimary = 1 AND so.zWarning = 0 AND so.class = 'GALAXY' -- gals w/ good spec
    AND (so.legacy_target1 & 32) > 0  -- LRGs

-- For QSOs:
INTO mydb.Spec_Photo_QSO
FROM SpecObj AS so
    JOIN PhotoObj AS po ON po.objID = so.bestObjID
    LEFT OUTER JOIN Photoz pz ON (pz.objID = so.bestObjID)
WHERE
    so.legacyPrimary = 1 AND so.zWarning = 0 AND so.class = 'QSO'  -- QSOs w/ good spec
```



The **deprecated** `gal` and `qso-old` tables were generated by this set of queries:

```sql
SELECT
    -- specObj metadata:
    so.specObjID, so.fluxObjId, so.bestObjID,
    so.sciencePrimary, so.sdssPrimary, so.legacyPrimary,
    so.firstRelease, so.survey,
    so.programname, so.mjd,
    so.plate, so.fiberID, so.run2d,
    so.ra, so.dec, so.cx, so.cy, so.cz,
    so.z, so.zErr,
    so.class, so.subClass, so.rChi2, so.DOF, so.rChi2Diff,
    so.z_person, so.class_person, so.comments_person,
    so.velDisp, so.velDispErr,
    so.waveMin, so.waveMax, so.wCoverage,
    so.spectroFlux_u, so.spectroFlux_g, so.spectroFlux_r, so.spectroFlux_i, so.spectroFlux_z,
    so.spectroFluxIvar_u, so.spectroFluxIvar_g, so.spectroFluxIvar_r, so.spectroFluxIvar_i, so.spectroFluxIvar_z,
    -- photometric catalog info:
    po.psfMag_r, po.psfMagErr_r, 
    -- modelMag uses best of deV or exp profiles (max like in r band):
    po.modelMag_u, modelMagErr_u,
    po.modelMag_g, modelMagErr_g, 
    po.modelMag_r, modelMagErr_r,
    po.modelMag_i, modelMagErr_i, 
    po.modelMag_z, modelMagErr_z,
    -- deV profile params:
    po.deVMag_u, po.deVMagErr_u, po.deVRad_u, po.deVAB_u,
    po.deVMag_g, po.deVMagErr_g, po.deVRad_g, po.deVAB_g,
    po.deVMag_r, po.deVMagErr_r, po.deVRad_r, po.deVAB_r,
    po.deVMag_i, po.deVMagErr_i, po.deVRad_i, po.deVAB_i,
    po.deVMag_z, po.deVMagErr_z, po.deVRad_z, po.deVAB_z,
    -- exp profile params:
    po.expMag_u, po.expMagErr_u, po.expRad_u, po.expAB_u,
    po.expMag_g, po.expMagErr_g, po.expRad_g, po.expAB_g,
    po.expMag_r, po.expMagErr_r, po.expRad_r, po.expAB_r,
    po.expMag_i, po.expMagErr_i, po.expRad_i, po.expAB_i,
    po.expMag_z, po.expMagErr_z, po.expRad_z, po.expAB_z,
    -- r-band model likelihoods:
    po.lnLDeV_r, po.lnLExp_r,
    -- petroMag can be used to estimate surface brightness
    po.petroMag_u, po.petroMagErr_u, po.petroRad_u, po.petroR50_u, po.petroR90_u,
    po.petroMag_g, po.petroMagErr_g, po.petroRad_g, po.petroR50_g, po.petroR90_g,
    po.petroMag_r, po.petroMagErr_r, po.petroRad_r, po.petroR50_r, po.petroR90_r,
    po.petroMag_i, po.petroMagErr_i, po.petroRad_i, po.petroR50_i, po.petroR90_i,
    po.petroMag_z, po.petroMagErr_z, po.petroRad_z, po.petroR50_z, po.petroR90_z,
    -- extinction corrections:
    po.extinction_u, po.extinction_g, po.extinction_r,
    po.extinction_i, po.extinction_z

-- For gal:
INTO mydb.Spec_Photo_Gal
FROM SpecObj AS so
    JOIN PhotoObj AS po ON so.bestObjID = po.objID
WHERE
    so.legacyPrimary = 1 AND so.zWarning = 0 AND so.class = 'GALAXY'

-- For qso-old:
INTO mydb.Spec_Photo_QSO_old
FROM SpecObj AS so
    JOIN PhotoObj AS po ON so.bestObjID = po.objID
WHERE
    so.legacyPrimary = 1 AND so.zWarning = 0 AND so.class = 'QSO'
```

```sql
SELECT
    -- specObj metadata:
    so.specObjID, so.fluxObjId, so.bestObjID,
    so.sciencePrimary, so.sdssPrimary, so.legacyPrimary,
    so.firstRelease, so.survey,
    so.programname, so.mjd,
    so.plate, so.fiberID, so.run2d,
    so.ra, so.dec, so.cx, so.cy, so.cz,
    so.z, so.zErr,
    so.class, so.subClass, so.rChi2, so.DOF, so.rChi2Diff,
    so.z_person, so.class_person, so.comments_person,
    so.velDisp, so.velDispErr,
    so.waveMin, so.waveMax, so.wCoverage,
    so.spectroFlux_u, so.spectroFlux_g, so.spectroFlux_r, so.spectroFlux_i, so.spectroFlux_z,
    so.spectroFluxIvar_u, so.spectroFluxIvar_g, so.spectroFluxIvar_r, so.spectroFluxIvar_i, so.spectroFluxIvar_z,
    -- photometric catalog info:
    po.psfMag_r, po.psfMagErr_r, 
    -- modelMag uses best of deV or exp profiles (max like in r band):
    po.modelMag_u, modelMagErr_u,
    po.modelMag_g, modelMagErr_g, 
    po.modelMag_r, modelMagErr_r,
    po.modelMag_i, modelMagErr_i, 
    po.modelMag_z, modelMagErr_z,
    -- deV profile params:
    po.deVMag_u, po.deVMagErr_u, po.deVRad_u, po.deVAB_u,
    po.deVMag_g, po.deVMagErr_g, po.deVRad_g, po.deVAB_g,
    po.deVMag_r, po.deVMagErr_r, po.deVRad_r, po.deVAB_r,
    po.deVMag_i, po.deVMagErr_i, po.deVRad_i, po.deVAB_i,
    po.deVMag_z, po.deVMagErr_z, po.deVRad_z, po.deVAB_z,
    -- exp profile params:
    po.expMag_u, po.expMagErr_u, po.expRad_u, po.expAB_u,
    po.expMag_g, po.expMagErr_g, po.expRad_g, po.expAB_g,
    po.expMag_r, po.expMagErr_r, po.expRad_r, po.expAB_r,
    po.expMag_i, po.expMagErr_i, po.expRad_i, po.expAB_i,
    po.expMag_z, po.expMagErr_z, po.expRad_z, po.expAB_z,
    -- r-band model likelihoods:
    po.lnLDeV_r, po.lnLExp_r,
    -- petroMag can be used to estimate surface brightness
    po.petroMag_u, po.petroMagErr_u, po.petroRad_u, po.petroR50_u, po.petroR90_u,
    po.petroMag_g, po.petroMagErr_g, po.petroRad_g, po.petroR50_g, po.petroR90_g,
    po.petroMag_r, po.petroMagErr_r, po.petroRad_r, po.petroR50_r, po.petroR90_r,
    po.petroMag_i, po.petroMagErr_i, po.petroRad_i, po.petroR50_i, po.petroR90_i,
    po.petroMag_z, po.petroMagErr_z, po.petroRad_z, po.petroR50_z, po.petroR90_z,
    -- extinction corrections:
    po.extinction_u, po.extinction_g, po.extinction_r,
    po.extinction_i, po.extinction_z
INTO mydb.Spec_Photo_QSO
FROM SpecObj AS so
    JOIN PhotoObj AS po ON so.bestObjID = po.objID
WHERE
    so.legacyPrimary = 1 AND so.zWarning = 0 AND so.class = 'QSO'
```



## License

This project is Copyright (c) 2020, 2021 by Tom Loredo and David Kent and licensed under the terms of the BSD 3-Clause license; see the included [LICENSE.txt](LICENSE.txt) file.


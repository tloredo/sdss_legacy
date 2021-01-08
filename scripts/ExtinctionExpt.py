"""
Experiment with extinction corrections for SEDs from the SDSS sample our
PZ group is using.

This uses a sample of SEDs provided by David Kent.

This script was developed using the following conda env:

    $ conda create --name py38astro python=3.8 anaconda
    $ pip install synphot
    $ pip install dust_extinction
    $ pip install dustmaps

Dust maps need to be fetched and stored in a suitable location; before running
this script; see the sed_sdss module docstring.

2021-01-03:  Refactored from a 2020-11 script using David Kent's data
"""
import matplotlib as mpl
from matplotlib import rc
from matplotlib.pyplot import *

from sdss_legacy import mgs, qso


rc('figure.subplot', left=0.1, bottom=.125, top=.95, right=.95)
# rc('lines', linewidth=2.0) # doesn't affect frame lines
rc('font', size=14)  # default for labels (not axis labels)
rc('font', family='serif')  # default for labels (not axis labels)
rc('axes', labelsize=18)
rc('xtick.major', pad=8)
rc('xtick', labelsize=14)
rc('ytick.major', pad=8)
rc('ytick', labelsize=14)
rc('savefig', dpi=300)
rc('axes.formatter', limits=(-4,4))


ion()


# Some objects to experiment with:
gals = [mgs.ith(i) for i in range(0, 700000, 70000)]
qsos = [qso.ith(i) for i in range(0, 70000, 7000)]


# Plot a set of SEDs on one figure, flux vs. log(wl).
fig = figure(figsize=(12,6))
xlabel(r'$\lambda$ (ang)')
ylabel(r'$F_\lambda$ ($10^{-17}$ erg/cm$^2$/s/ang)')
for obj in gals:
    semilogx(obj.sed.lam, obj.sed.flux, '-')
axhline(0, c='k', lw=1)
wticks = [4000, 5000, 6000, 7000, 8000, 9000]
xticks(wticks)
gca().get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())


# Plot some adjusted SEDs, for different choices of map and ext. curves.

def plot_dered(axnum, sed, curves, map, Rv=3.1):
    """
    Plot a dereddened spectrum as a curve, and for each point, an error
    bar dipping to the original (unadjusted) spectrum value.
    """

    if axnum in (3,4):
        xlabel(r'$\lambda$ (ang)')
    if axnum in (1,3):
        ylabel(r'$F_\lambda$ ($10^{-17}$ erg/cm$^2$/s/ang)')

    aflux, aivar, fac = sed.dereddened(curves=curves, map=map, factor=True)

    # Using a fill between adjusted and unadjusted spectra looks wierd!
    # fill_between(sed.lam, sed.flux, aflux)
    # xscale('log')

    # Plot the adjusted flux, and error bars down to the raw flux.
    semilogx(sed.lam, aflux, 'k-', lw=1, alpha=.5)
    mid = 0.5*(sed.flux + aflux)
    hdeltas = 0.5*(aflux - sed.flux)
    errorbar(sed.lam, mid, yerr=hdeltas, fmt='none', ecolor='C0', alpha=.5)
    # axhline(0, c='k', lw=1)
    y_l, y_u = ylim()
    ylim(0, y_u)

    # Plot extinction factor against rt axis.
    ax = gca()
    ax.tick_params(direction='in')
    xticks(wticks)
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

    rax = ax.twinx()
    rax.tick_params(direction='in')
    rax.plot(sed.lam, fac, '-', color='firebrick', lw=2)
    if axnum in (2,4):
        rax.set_ylabel('Extinction factor', color='firebrick')
    rax.set_ylim(0., 1.)

    annot = '{}, {}'.format(curves, map)
    text(0.05, 0.75, annot, horizontalalignment='left',
         verticalalignment='bottom', transform=ax.transAxes)


fig = figure(figsize=(15,8))
fig.subplots_adjust(right=0.9)  # for rt axis=
sed = gals[0].sed
subplot(221)
plot_dered(1, sed, 'CCM89', 'SFD')
subplot(222)
plot_dered(2, sed, 'CCM89', 'Planck')
subplot(223)
plot_dered(3, sed, 'F99', 'SFD')
subplot(224)
plot_dered(4, sed, 'F99', 'Planck')

fig = figure(figsize=(15,8))
fig.subplots_adjust(right=0.9)  # for rt axis=
sed = gals[7].sed
subplot(221)
plot_dered(1, sed, 'CCM89', 'SFD')
subplot(222)
plot_dered(2, sed, 'CCM89', 'Planck')
subplot(223)
plot_dered(3, sed, 'F99', 'SFD')
subplot(224)
plot_dered(4, sed, 'F99', 'Planck')

fig = figure(figsize=(15,8))
fig.subplots_adjust(right=0.9)  # for rt axis=
sed = qsos[3].sed
subplot(221)
plot_dered(1, sed, 'CCM89', 'SFD')
subplot(222)
plot_dered(2, sed, 'CCM89', 'Planck')
subplot(223)
plot_dered(3, sed, 'F99', 'SFD')
subplot(224)
plot_dered(4, sed, 'F99', 'Planck')

galex_motion
============

`galex_motion` is a module to make grabbing and coadding measurements of GALEX fluxes and upper limits from the master source catalog (MCAT) en masse easy. The code will propagate the proper motion  of the source(s) you want to retrieve fluxes for to the appropriate epoch of each GALEX exposure that might contain the source. Hence the "motion" in the name. It also combines FUV and NUV measurements that are sometimes erroneously identified as separate FUV and NUV sources in the MCAT. For targets with matches in multiple GALEX exposures, outliers are rejected to mitigate the influence of flares. There is also a utility to correct fluxes in the detectors nonlinear regime.

For a module with much greater power and flexibility, including the ability to create intra-exposure lightcurves and access the photon-level GALEX data, use gPhoton instead (https://github.com/cmillion/gPhoton). Alternatively, if you just want a quick look at UV fluxes for a single target, the GALEXView interface is likely the friendliest option (http://galex.stsci.edu/galexView/).

If this is your first time playing with GALEX data, have a look at Morrisey+ 2007 (https://ui.adsabs.harvard.edu/abs/2007ApJS..173..682M).

This code uses the `requests`, `numpy`, and `astropy` modules.

# Quick Start

Let's get the GALEX fluxes for Epsilon Eridani, a bright, active K star with a high  proper motion.

```
import galex_motion as gm

# Eps Eri's coordinates and proper motions from SIMBAD
ra, dec = 53.2326873475, -9.4582586569
pm_ra, pm_dec = -975.17, 19.49

nuv, fuv = gm.extract_and_coadd(ra, dec, pm_ra, pm_dec)
# nuv and fuv are (flux, error) tuples
# fluxes are in counts s-1 -- see Morrisey+ 2007 for conversion to magnitude

print(nuv[0])
# 392.68309999999997

# the nuv flux `nuv[0]`, is > 300, which is in the range where the GALEX
# response becomes nonlinear. Better correct it.
corrected_nuv = gm.nonlinearity_correction(*nuv, band='nuv')

print(corrected_nuv[0])
# 504.5195387279937
```

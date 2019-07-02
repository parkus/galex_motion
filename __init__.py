# the SQL query for this code was modified from one in the gPhoton package
# (https://github.com/cmillion/gPhoton).
#
# description of MCAT columns available at
# http://www.galex.caltech.edu/wiki/Public:Documentation/Appendix_A.1


import numpy as np
import requests
from astropy import table, time, units as u
_j2000_jd = 2451545.0
_usable_fov = 1.2
_aper_7_radius = 23 / 2.  # pix
_aper_7_area = np.pi * (_aper_7_radius * 1.5) ** 2  # arcsec2

baseURL = ('https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/Galex'
           'PhotonListQueryTest?query=')
MCATDB = 'GR6Plus7.dbo'
default_columns = "ra, dec, NUV_FLUX_APER_7, NUV_FLUXERR_APER_7, " \
                  "FUV_FLUX_APER_7, FUV_FLUXERR_APER_7, " \
                  "nuv_skybg, fuv_skybg, vpe.fexptime, " \
                  "vpe.nexptime, nuv_artifact, fuv_artifact, " \
                  "vpe.fexpstar,  vpe.fexpend, vpe.nexpstar, vpe.nexpend, " \
                  "fov_radius".split(', ')


def _combine_columns(addnl_columns):
    cols = default_columns
    if addnl_columns is not None:
        cols += addnl_columns
    return cols


def _assemble_http_query(ra, dec, radius, addnl_columns=None):
    cols = _combine_columns(addnl_columns)
    cols_str = ', '.join(cols)
    query = "{baseURL}select " \
            "{cols} " \
            "from {MCATDB}.visitphotoobjall as vpo " \
            "" \
            "inner join {MCATDB}.visitphotoextract as vpe on " \
            "vpo.photoextractid=vpe.photoextractid " \
            "" \
            "inner join {MCATDB}.fGetNearbyVisitObjEq({ra},{dec},{radius}) " \
            "as nb on vpo.objid=nb.objid " \
            "" \
            "&format=extjs" \
            "".format(baseURL=baseURL, MCATDB=MCATDB, ra=ra, dec=dec,
                      radius=radius*60., cols=cols_str)
    return query


def _strip_prefix(name):
    pieces = name.split('.')
    return pieces[-1]


def compute_expected_positions(ra, dec, pm_ra, pm_dec, jd):
    pm = np.array((pm_ra, pm_dec))/3600/1000./365 # deg/d
    def translate(x, dxdt):
        return x + dxdt*(jd-_j2000_jd)
    return map(translate, (ra, dec), pm)


def _dist(xy0, xy1):
    xy0, xy1 = [np.reshape(a, (-1, 2)) for a in (xy0, xy1)]
    return np.sqrt(np.sum((xy1 - xy0)**2, axis=1))


def _mcattime2jd(time_mcat):
    time_unix = time.Time(time_mcat, format='unix')
    return time_unix.jd





def fetch_mcat_data(ra, dec, search_radius=2 / 60., addnl_columns=None,
                    timeout=60., connect_retries=10):

    # for reference, 1-sigma astrometric uncertainty is 0.4"
    query = _assemble_http_query(ra, dec, search_radius,
                                 addnl_columns=addnl_columns)
    count = 1
    while count <= 10:
        try:
            r = requests.get(query, timeout=timeout)
            break
        except requests.ConnectionError:
            count += 1
            continue

    try:
        rows = r.json()['data']['Tables'][0]['Rows']
    except:
        raise ValueError('Hmmm no data returned, instead got:\n{}'
                         ''.format(r.content))

    cols = _combine_columns(addnl_columns)
    names = map(_strip_prefix, cols)
    if len(rows) == 0:
        tbl = table.Table(names=names, masked=True)
    else:
        tbl = table.Table(rows=rows, names=names, masked=True)
    tbl.meta['search_position'] = ra, dec
    tbl.meta['search_radius'] = search_radius
    return tbl


def add_position_offset_column(ra, dec, pm_ra, pm_dec, tbl):
    # I suppose computing the exposure midpoint would be best, but the target
    # motion should always be negligible relative to exposure time so it
    # isn't necessary
    fuv_and_nuv_times = np.array((tbl['nexpstar'], tbl['fexpstar']))
    times = np.max(fuv_and_nuv_times, axis=0)
    jd = _mcattime2jd(times)
    reference_position = compute_expected_positions(ra, dec, pm_ra, pm_dec, jd)
    reference_position = np.array(reference_position).T
    source_positions = np.array((tbl['ra'], tbl['dec'])).T
    tbl['offset'] = _dist(reference_position, source_positions)


def extract_source(ra, dec, pm_ra, pm_dec, match_radius=4./3600.,
                   search_radius=25./60, query_timeout=60.,
                   upper_limits=True):
    # compute coordinates at start and end of galex lifetime (2004-2011)
    def correct_coords(year):
        jd = (year - 2000)*365.25 + _j2000_jd
        return compute_expected_positions(ra, dec, pm_ra, pm_dec, jd)
    ra04, dec04 = correct_coords(2004)
    ra08, dec08 = correct_coords(2008)
    ra12, dec12 = correct_coords(2012)

    # enlarge radius to catch source over that range of positions
    path_length = np.sqrt((dec12 - dec04)**2 + (ra12 - ra04)**2)
    capture_radius = path_length/2 + search_radius

    # pull MCAT data from larger area to catch any tiles that might have had
    # the source
    tbl = fetch_mcat_data(ra08, dec08, capture_radius, timeout=query_timeout)

    # add column of offsets from expected source positions
    add_position_offset_column(ra, dec, pm_ra, pm_dec, tbl)

    # match source and get fluxes for each band
    def get_fluxes(band):
        return get_nearest_source_fluxes(tbl, band, match_radius=match_radius,
                                         upper_limits=upper_limits)

    return get_fluxes('nuv'), get_fluxes('fuv')


def extract_and_coadd(ra, dec, pm_ra, pm_dec, match_radius=4./3600.,
                      search_radius=25./60, sigma_clip=3., query_timeout=60.,
                      upper_limits=True):
    data = extract_source(ra, dec, pm_ra, pm_dec, match_radius,
                          search_radius, query_timeout, upper_limits)
    nuv_data, fuv_data = data
    return coadd_fluxes(*nuv_data[:3], sigma_clip=sigma_clip), \
           coadd_fluxes(*fuv_data[:3], sigma_clip=sigma_clip)


def coadd_fluxes(fluxes, errs, expts, sigma_clip=3.):
    if len(fluxes) == 0:
        return -999., -999.
    fluxes, errs, expts = map(np.asarray, (fluxes, errs, expts))
    if np.all(fluxes <= -99.):
        # just pick most restrictive limit
        return -999., np.min(errs)
    else:
        def apply_filter(keep):
            return [a[keep] for a in (fluxes, errs, expts)]
        keep = fluxes > 0
        fluxes, errs, expts = apply_filter(keep)
        median = np.median(fluxes)
        keep = (fluxes - median)/errs < sigma_clip
        fluxes, errs, expts = apply_filter(keep)
        coadd_flux = np.sum(fluxes*expts)/np.sum(expts)
        coadd_err = np.sqrt(np.sum((errs*expts)**2))/np.sum(expts)
        return coadd_flux, coadd_err


def get_nearest_source_fluxes(tbl, band, match_radius=2/3600.,
                              upper_limits=True):
    if len(tbl) == 0:
        return [], [], [], []

    # for each unique time (i.e. each unique exposure), find nearest source
    # note that there  can be multiple exposures  for a single tile, so don't
    # use tile id for this
    letter = band[0]
    fluxcol = band.upper() + '_FLUX_APER_7'
    errcol = band.upper() + '_FLUXERR_APER_7'

    unique_expend = np.unique(tbl[letter + 'expend'])
    unique_expend = unique_expend[unique_expend > 0]
    fluxes, errs, expts, offsets = [], [], [], []
    for expend in unique_expend:
        exp_tbl = tbl[tbl[letter + 'expend'] == expend]
        i_closest = np.argmin(exp_tbl['offset'])
        flux = exp_tbl[fluxcol][i_closest]
        err = exp_tbl[errcol][i_closest]

        too_far = exp_tbl['offset'][i_closest] > match_radius
        null_value = flux <= -99.
        if too_far or null_value:
            # first make sure tile center is close enough that source would have
            # been within the match radius if it was present. there is not a column
            # in the mcat database for boresight position, but if the center is
            # close enough then the distance between the points and the search
            # coordinates plust the points and the boresight will all be greater
            # than the detector radius
            source_pts = np.array((exp_tbl['ra'], exp_tbl['dec'])).T
            search_pt = np.array(tbl.meta['search_position'])
            dist = _dist(source_pts, search_pt)
            if np.all(dist + exp_tbl['fov_radius'] > _usable_fov/2.):
                continue

            if upper_limits:
                expt = np.median(exp_tbl[letter + 'exptime'])
                bg_cps = np.mean(exp_tbl[band + '_skybg']) * _aper_7_area
                bg_cpserr = np.sqrt(bg_cps/expt)
                lim = 2*bg_cpserr
                fluxes.append(-999.)
                errs.append(lim)
                expts.append(expt)
            else:
                fluxes.append(-999.)
                errs.append(-999.)
                expts.append(-999.)
            offsets.append(-999.)
        else:
            fluxes.append(flux)
            errs.append(err)
            offsets.append(exp_tbl['offset'][i_closest])
            expts.append(exp_tbl[letter + 'exptime'][i_closest])

    return fluxes, errs, expts, offsets


def nonlinearity_correction(cps, err, band):
    """Good rule of thumb is not to trust correction when error following
    correction is >3x the error prior to correction."""

    # I fit these from the plot for the APER_7 curve to match the catalog
    # values this code pulls
    if band.lower() == "fuv":
        a, b, c0 = -0.23428377,  1.69610651, -0.51978853
    elif band.lower() == "nuv":
        a, b, c0 = -0.15198511,  1.60029826, -0.621034
    else:
        raise ValueError("Band must be NUV or FUV.")

    # find where the corrected values are closest to the measured values
    logPR0 = np.sqrt(c0/a)
    logMR0 = a*logPR0**2 + b*logPR0 + c0

    logMR = np.log10(cps)
    c = c0 - logMR

    logPR = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
    PR = 10**logPR
    dlogPR_dc = 1./np.sqrt(b**2 - 4*a*c)
    dc_dPR = 1/np.log(10)/cps
    PRerr = 10**logPR*np.log(10)*dlogPR_dc*dc_dPR*err

    correct = logMR > logMR0
    PR[~correct] = cps[~correct]
    PRerr[~correct] = err[~correct]
    return PR, PRerr



# the SQL query for this code was modified from one in the gPhoton package
# (https://github.com/cmillion/gPhoton).
#
# description of MCAT columns available at
# http://www.galex.caltech.edu/wiki/Public:Documentation/Appendix_A.1


import numpy as np
import requests
from astropy import table, time
_j2000_jd = 2451545.0
_usable_fov = 1.2


baseURL = ('https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/Galex'
           'PhotonListQueryTest?query=')
MCATDB = 'GR6Plus7.dbo'
default_columns = "ra, dec, nuv_flux, nuv_fluxerr, fuv_flux, fuv_fluxerr, " \
                  "vpe.fexptime, vpe.nexptime, nuv_artifact, fuv_artifact, " \
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


def circle_intersections(xyr0, xyr1):
    """
    Compute the points where two circles interesect.

    Parameters
    ----------
    circles : 2x3 array-like
        [[x0,y0,r0],[x1,y1,r1]] - the radius and center coordinates of the two
        circles

    Returns
    -------
    xpts : list or array
        If the circles do not interset, an empty list is returned. If one
        circle is enclosed in another, the index ([0] or [1]) of that circle is
        returned. Otherwise, the intersection points are returned as the array
        [[xi0,yi0], [xi1,yi1]] with xi0 >= xi1. In the rare case that the
        just touch, the one intersection point will be returned twice.
    """
    (x0, y0, r0), (x1, y1, r1) = xyr0, xyr1
    d = _dist((x0, y0), (x1, y1))
    if np.isclose(d, (r0 + r1)):
        frac = r0/r1
        xi = x0 + frac*(x1 - x0)
        yi = y0 + frac*(y1 - y0)
        return (xi, yi),
    if d > (r0 + r1):  # if the circles do not intersect
        return ()
    elif d < abs(r1 - r0):  # if one circle is within another
        return ()
    else:
        # compute intersection. some variables are arbitrarly assigned just to
        # make the math more concise. the math may be worked out by
        # simultaneously solving the equations for two circles
        q = (r0 ** 2 - r1 ** 2 + x1 ** 2 - x0 ** 2 + y1 ** 2 - y0 ** 2) / 2.0
        dx, dy = (x1 - x0), (y1 - y0)
        a = 1 + dx ** 2 / dy ** 2
        b = -2 * x0 - 2 * q * dx / dy ** 2 + 2 * y0 * dx / dy
        c = x0 ** 2 + y0 ** 2 + q ** 2 / dy ** 2 - 2 * q * y0 / dy - r0 ** 2
        xi0 = (-b + np.sqrt(b ** 2 - 4 * a * c)) / 2 / a
        xi1 = (-b - np.sqrt(b ** 2 - 4 * a * c)) / 2 / a
        yi0, yi1 = (q - xi0 * dx) / dy, (q - xi1 * dx) / dy
        return (xi0, yi0), (xi1, yi1)


def _infer_center(xyr0, xyr1, xyr2):
    ixs = circle_intersections(xyr0, xyr1)

    if len(ixs) == 0:
        raise ValueError("Uh oh, there doesn't seem to be a common center.")
    elif len(ixs)  == 1:
        # if somehow two points are exactly opposed, the circles defined by
        # their radius from the center will just graze and give one point
        return ixs[0]
    else:
        # otherwise, there will be two points of intersection. find the
        # actual center by using a point to find which is closer to the
        # circle its radius defines
        xs, ys = np.array(ixs).T
        x, y, r = xyr2
        radii = np.sqrt((xs - x) ** 2 + (ys - y) ** 2)
        error = np.abs(radii - r)
        i = np.argmin(error)
        if error[i] > 0.001:
            raise ValueError("Hmmmm, that ain't right.")
        return xs[i], ys[i]


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
                    timeout=60.):

    # for reference, 1-sigma astrometric uncertainty is 0.4"
    query = _assemble_http_query(ra, dec, search_radius,
                                 addnl_columns=addnl_columns)
    r = requests.get(query, timeout=timeout)
    try:
        rows = r.json()['data']['Tables'][0]['Rows']
    except:
        raise ValueError('Hmmm no data returned, instead got:\n{}'
                         ''.format(r.content))

    cols = _combine_columns(addnl_columns)
    names = map(_strip_prefix, cols)
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


def match_source(ra, dec, pm_ra, pm_dec, match_radius=4./3600.,
                 search_radius=25./60, sigma_clip=3., query_timeout=60.,
                 include_coarse_limits=True):
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
                                         include_coarse_limits=include_coarse_limits)
    fuv_flux, fuv_err = coadd_fluxes(*get_fluxes('fuv'), sigma_clip=sigma_clip)
    nuv_flux, nuv_err = coadd_fluxes(*get_fluxes('nuv'), sigma_clip=sigma_clip)

    return fuv_flux, fuv_err, nuv_flux, nuv_err


def match_and_coadd(ra, dec, pm_ra, pm_dec, match_radius=4./3600.,
                    search_radius=25./60, sigma_clip=3., query_timeout=60.,
                    include_coarse_limits=True)):
    pass


def coadd_fluxes(fluxes, errs, expts, sigma_clip=3.):
    fluxes, errs, expts = map(np.asarray, (fluxes, errs, expts))
    if np.all(np.isnan(fluxes)):
        # just pick most restrictive limit
        return np.nan, np.min(errs)
    else:
        def apply_filter(keep):
            return [a[keep] for a in (fluxes, errs, expts)]
        keep = ~np.isnan(fluxes)
        fluxes, errs, expts = apply_filter(keep)
        median = np.median(fluxes)
        keep = (fluxes - median)/errs < sigma_clip
        fluxes, errs, expts = apply_filter(keep)
        coadd_flux = np.sum(fluxes*expts)/np.sum(expts)
        coadd_err = np.sqrt(np.sum((errs*expts)**2))/np.sum(expts)
        return coadd_flux, coadd_err


def get_nearest_source_fluxes(tbl, band, match_radius=2/3600.,
                              include_coarse_limits=True):
    # for each unique time (i.e. each unique exposure), find nearest source
    # note that there  can be multiple exposures  for a single tile, so don't
    # use tile id for this
    letter = band[0]

    unique_expend = np.unique(tbl[letter + 'expend'])
    unique_expend = unique_expend[unique_expend > 0]
    fluxes, errs, expts = [], [], []
    for expend in unique_expend:
        exp_tbl = tbl[tbl[letter + 'expend'] == expend]
        if not np.any(exp_tbl['offset'] < match_radius):
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
            if include_coarse_limits and len(exp_tbl) > 10:
                other_fluxes = exp_tbl[letter + 'uv_flux']
                other_errs = exp_tbl[letter + 'uv_fluxerr']
                a, b = np.polyfit(np.log(other_fluxes), np.log(other_errs), 1)

                # twice the error for which  S/N = 1 (log(err) = log(flux))
                log_err_SN1 = b/(1-a)
                lim = 2*np.exp(log_err_SN1)

                fluxes.append(np.nan)
                errs.append(lim)
            else:
                fluxes.append(np.nan)
                errs.append(np.nan)
        else:
            i_closest = np.argmin(exp_tbl['offset'])
            flux = exp_tbl[letter + 'uv_flux'][i_closest]
            err = exp_tbl[letter + 'uv_fluxerr'][i_closest]
            fluxes.append(flux)
            errs.append(err)
        expts.append(exp_tbl[letter + 'exptime'][i_closest])

    return fluxes, errs, expts

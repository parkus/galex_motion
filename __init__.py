# the SQL query for this code was modified from one in the gPhoton package
# (https://github.com/cmillion/gPhoton).
#
# description of MCAT columns available at
# http://www.galex.caltech.edu/wiki/Public:Documentation/Appendix_A.1


import numpy as np
import requests
from astropy import table, time, units as u


baseURL = ('https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/Galex'
           'PhotonListQueryTest?query=')
MCATDB = 'GR6Plus7.dbo'
default_columns = "ra, dec, nuv_mag, fuv_mag, vpe.fexptime, " \
                  "vpe.nexptime, nuv_artifact, fuv_artifact, vpe.fexpstar, " \
                  "vpe.fexpend, vpe.nexpstar, vpe.nexpend".split(', ')


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
                      radius=radius, cols=cols_str)
    return query


def _strip_prefix(name):
    pieces = name.split('.')
    return pieces[-1]


def _move(ra, dec, pm_ra, pm_dec, year):
    translate = lambda x, dxdt: x + dxdt*(year-2000)
    return map(translate, (ra, dec), (pm_ra, pm_dec))


def _dist(xy0, xy1):
    xy0, xy1 = [np.reshape(a, (-1,2)) for a in (xy0, xy1)]
    return np.sqrt(np.sum(xy1**2 - xy0**2, axis=1))


def _galextime2year(time_mcat):
    time_unix = time.Time(time_mcat, format='unix')
    time_j2000 = time.Time(2451545.0, format='jd')
    time_year = time_unix - time_j2000
    return time_year.to(u.yr).value


def fetch_mcat_data(ra, dec, pm_ra=0, pm_dec=0, search_radius=2 / 60.,
                  addnl_columns=None, timeout=60.):

    # region no source motion
    if pm_ra == 0 and pm_dec == 0:
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
        coords = np.array((ra, dec))
        tbl['offset'] = _dist((tbl['ra'], tbl['dec']), coords[None,:])
        return tbl
    # endregion

    # region source is moving
    ## compute coordinates at start and end of galex lifetime (2004-2011)
    correct_coords = lambda year: _move(ra, dec, pm_ra, pm_dec, year)
    ra04, dec04 = correct_coords(2004)
    ra08, dec08 = correct_coords(2008)
    ra12, dec12 = correct_coords(2012)

    ## enlarge radius to catch source over that range of positions
    path_length = np.sqrt((dec12 - dec04)**2 + (ra12 - ra04)**2)
    capture_radius = (path_length*60)/2 + search_radius

    ## search over the larger area to catch any  tiles that might have had the
    ## source
    tbl = fetch_mcat_data(ra08, dec08, 0, 0, capture_radius, timeout=timeout)

    ## compute coordinates at each time and re-reference the offset column
    ## I suppose computing the exposure midpoint would be best, but the target
    ## motion should always be negligible relative to exposure time so it
    ## isn't necessary
    fuv_and_nuv_times = np.array((tbl['nexpstar', 'fexpstar']))
    times = np.max(fuv_and_nuv_times, axis=0)
    years = _galextime2year(times)
    corrected_coords = correct_coords(ra, dec, pm_ra, pm_dec, years)
    tbl['offset'] = _dist(tbl['ra'], tbl['dec'], corrected_coords)

    return  tbl


def coadd_mags(rad, dec, pm_ra=0, pm_dec=0, search_radius=2/60.,
               backout_limit=True, mask_flags=0b1111111111):
    for band in ['fuv', 'nuv']:
        fuv_flags = tbl[band + '_artifact']
        mask = np.bitwise_and(fuv_flags, mask_flags)
        tbl[band + '_mag'].mask = mask


def bootstrap_limit():
    pass


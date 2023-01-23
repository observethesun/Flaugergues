"""Coordinate reconstruction tools."""
from functools import partial
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from astropy.coordinates import cartesian_to_spherical
from astropy.time import Time
from sunpy.coordinates.sun import L0, B0, orientation, P
import ephem

def _from_local(row):
    """Local time to utc."""
    rec_t = row.recorded_time
    if row.rising:
        if rec_t.hour < 12:
            delta = rec_t - (pd.to_datetime(rec_t.date()) + pd.Timedelta(hours=12))
        else:
            delta = rec_t - (pd.to_datetime(rec_t.date()) + pd.Timedelta(hours=24))
    else:
        if rec_t.hour < 12:
            delta = rec_t - pd.to_datetime(rec_t.date())
        else:
            delta = rec_t - (pd.to_datetime(rec_t.date()) + pd.Timedelta(hours=12))
    return row.utc_noon + delta

def _from_sidereal(row):
    """Sidereal time to utc."""
    sid2norm = 1 / 1.0027378836477503
    dt = (row.recorded_time - row.sidereal_at_noon) * sid2norm + row.utc_noon
    time = dt.time()
    return (pd.to_datetime(row.recorded_time.date()) +
            pd.Timedelta(hours=time.hour,
                         minutes=time.minute,
                         seconds=time.second,
                         microseconds=time.microsecond))

def _get_time(row):
    """Mapping selector."""
    return _from_sidereal(row) if row.sidereal else _from_local(row)

def rectime2utctime(rec, obs_loc):
    """Map recorded_time to utc_time."""
    obs = ephem.Observer()
    obs.lat, obs.long = str(obs_loc.lat.deg), str(obs_loc.lon.deg)
    sun = ephem.Sun()
    rec['utc_noon'] = (rec.recorded_time
                       .apply(lambda x: obs.next_transit(sun, start=pd.to_datetime(x.date())).datetime()))
    t = Time(rec.utc_noon, scale='utc', location=obs_loc)
    rec['sidereal_at_noon'] = t.sidereal_time('apparent')
    rec.sidereal_at_noon = (rec.sidereal_at_noon.apply(lambda x: pd.Timedelta(hours=x)) +
                            pd.to_datetime(rec.recorded_time.dt.date))
    rec['utc_time'] = rec.apply(_get_time, axis=1)
    rec = rec.drop(columns=['utc_noon', 'sidereal_at_noon'])
    return rec

def rot2d(vec, theta, deg=False):
    """Rotate 2D vector by angle theta counter-clockwise."""
    vec = np.asarray(vec)
    shape = vec.shape
    theta = np.deg2rad(theta) if deg else theta
    c, s = np.cos(theta), np.sin(theta)
    return (np.array([[c, -s], [s, c]]) @ vec.reshape(2, 1)).reshape(shape)

def xy_to_hc(vec, B_0, L_0, P_0, deg=True):
    """Map x, y coordinates on the unit solar disk from the coordinate system aligned
    with the motion of the Sun to heliographic coordinates."""
    vec = np.asarray(vec)
    x, y = rot2d(vec, -P_0, deg)
    vec3 = np.array([np.sqrt(1 - x*x - y*y), x, y])
    vec3 = Rotation.from_euler('y', -B_0, deg).apply(vec3)
    _, lat, long = cartesian_to_spherical(*vec3)

    if deg:
        lat, car_long = lat.deg, (long.deg + L_0) % 360
        cmd = long.deg - 360 if long.deg > 90 else long.deg
    else:
        lat, car_long = lat.rad, (long.rad + L_0) % 2*np.pi
        cmd = long.rad - 2*np.pi if long.deg > np.pi else long.rad

    return {'lat': lat,
            'long': car_long,
            'cmd': cmd,
            'x_hgc': vec3[1],
            'y_hgc': vec3[2]}

def get_coords_ob(rec, alpha=45, pole=None, deg=True):
    """Compute x, y coordinates from observations with oblique wires.
    The x-axis is aligned with the line of motion of the Sun, the y-axis points
    to the northern hemisphere, and the radius of the Sun is one.

    Parameters
    ----------
    rec : pandas.DataFrame
        Record with contact times.
    alpha : float
        Inclination of the oblique wires.
    pole : str, optional
        Indicates the pole that contacts the horizontal wire. 'S' or 'IN' for the south pole,
        'N' or 'SP' for the north pole. If not given, the value from the record is used.
    deg : bool
        If True, angles are degrees, else radians. Default True.

    Returns
    -------
    out : list
        Dictionary with x, y coordinates and additional parameters for each sunspot in the record.
    """
    if set(rec.event.unique()) != {'hr', 'ob'}:
        raise ValueError('Invalid record {}.'.format(rec.obs_id.iloc[0]))
    spot_names = list(rec.object.loc[rec.object != 'S'].unique())
    sun = rec.loc[(rec.object == 'S') & (rec.event == 'hr')]
    if len(sun) != 2:
        raise ValueError('Expected 2 contact times for the Sun, found {} in record {}.'
                         .format(len(sun), rec.obs_id.iloc[0]))

    if pole is None:
        pole = sun.pole.iloc[0]
    sign = 1 if pole in ['S', 'IN'] else -1
    complemented = False
    err = 0
    pi = np.pi / 2
    if deg:
        alpha = np.deg2rad(alpha)

    out = []
    d = (sun.utc_time.iloc[1] - sun.utc_time.iloc[0]).seconds
    for name in spot_names:
        spot = rec.loc[rec.object == name]
        if set(spot.event.values) != {'hr', 'ob'}:
            print('Incomplete spot {} in record {}.'.format(name, rec.obs_id.iloc[0]))
            continue
        p1 = (spot.loc[spot.event == 'hr', 'utc_time'].iloc[0] - sun.utc_time.iloc[0]).seconds
        if len(spot.loc[spot.event == 'ob']) == 2:
            d1 = (spot.loc[spot.event == 'hr', 'utc_time'].iloc[0] -
                  spot.loc[spot.event == 'ob', 'utc_time'].iloc[0]).seconds
            d2 = (spot.loc[spot.event == 'ob', 'utc_time'].iloc[1] -
                  spot.loc[spot.event == 'hr', 'utc_time'].iloc[0]).seconds
        else:
            complemented = True
            if (spot.loc[spot.event == 'hr', 'utc_time'].iloc[0] >
                spot.loc[spot.event == 'ob', 'utc_time'].iloc[0]):
                d1 = (spot.loc[spot.event == 'hr', 'utc_time'].iloc[0] -
                      spot.loc[spot.event == 'ob', 'utc_time'].iloc[0]).seconds
                d2 = d1
            else:
                d2 = (spot.loc[spot.event == 'ob', 'utc_time'].iloc[0] -
                      spot.loc[spot.event == 'hr', 'utc_time'].iloc[0]).seconds
                d1 = d2
        p2 = p1 - d1
        if d1 != d2:
            err = np.arctan(np.tan(alpha) * (d2 - d1) / (d1 + d2))
            r_sun = d * np.cos(err) / 2
            ratio = np.tan(alpha + err) / np.tan(pi - err)
            x = r_sun / np.tan((pi + err) / 2) - (p1 + p2 * ratio) / (1 + ratio)
            y = -sign * np.tan(alpha + err) * (x - r_sun / np.tan((pi + err) / 2) + p2)
        else:
            r_sun = d / 2
            x = r_sun - p1
            y = sign * np.tan(alpha) * (p1 - p2)

        y = y / r_sun - sign
        x /= r_sun

        out.append({'name': name,
                    'x': x,
                    'y': y,
                    'r_sun': r_sun,
                    'err': np.rad2deg(err) if deg else err,
                    'complemented': complemented})
    return out

def get_coords_hv(rec, incl=None, deg=True):
    """Compute x, y coordinates from observations with horizontal and vertical wires.
    The x-axis is aligned with the line of motion of the Sun, the y-axis points
    to the northern hemisphere, and the radius of the Sun is one.

    Parameters
    ----------
    rec : pandas.DataFrame
        Record with contact times.
    incl : float, optional
        Inclination to be used for incomplete records.
    deg : bool
        If True, angles are degrees, else radians. Default True.

    Returns
    -------
    out : list
        Dictionary with x, y coordinates and additional parameters for each sunspot in the record.
    """
    out = []
    if set(rec.event.unique()) != {'h', 'v'}:
        raise ValueError('Invalid record {}.'.format(rec.obs_id.iloc[0]))
    spot_names = list(rec.object.loc[rec.object != 'S'].unique())
    sh = rec.loc[(rec.object == 'S') & (rec.event == 'h')]
    sv = rec.loc[(rec.object == 'S') & (rec.event == 'v')]
    rising = rec.rising.iloc[0]
    complemented = False

    if len(sh) == 0 or len(sv) == 0:
        raise ValueError('Incomplete Sun in record {}.'.format(rec.obs_id.iloc[0]))
    if len(sh) != 2 and len(sv) != 2:
        raise ValueError('Incomplete Sun in record {}.'.format(rec.obs_id.iloc[0]))

    if len(sh) == 2 and len(sv) == 2:
        dh = (sh.utc_time.iloc[1] - sh.utc_time.iloc[0]).seconds
        dv = (sv.utc_time.iloc[1] - sv.utc_time.iloc[0]).seconds
        incl = np.arctan(dv / dh)
        if rising:
            incl = -incl
        r_sun = dh * np.sin(np.abs(incl)) / 2
    else:
        if incl is None:
            raise ValueError('Missing incl to process the incomplete record {}.'
                             .format(rec.obs_id.iloc[0]))
        complemented = True
        incl = np.deg2rad(incl) if deg else incl
        if len(sh) == 2:
            dh = (sh.utc_time.iloc[1] - sh.utc_time.iloc[0]).seconds
            r_sun = dh * np.sin(np.abs(incl)) / 2
            dv = dh * np.tan(np.abs(incl))
        if len(sv) == 2:
            dv = (sv.utc_time.iloc[1] - sv.utc_time.iloc[0]).seconds
            dh = dv / np.tan(np.abs(incl))
            r_sun = dh * np.sin(np.abs(incl)) / 2

    for name in spot_names:
        spot = rec.loc[rec.object == name]
        if len(spot) != 2:
            print('Incomplete spot {} in record {}.'.format(name, rec.obs_id.iloc[0]))
            continue
        if set(spot.event.values) != {'h', 'v'}:
            print('Incomplete spot {} in record {}.'.format(name, rec.obs_id.iloc[0]))
            continue
        if spot.loc[spot.event == 'v', 'utc_time'].iloc[0] > sv.utc_time.iloc[0]:
            d = (spot.loc[spot.event == 'v', 'utc_time'].iloc[0] - sv.utc_time.iloc[0]).seconds
            x = 1 - 2*d / dv
        else:
            d = (sv.utc_time.iloc[0] - spot.loc[spot.event == 'v', 'utc_time'].iloc[0]).seconds
            x = -1 + 2*d / dv
        if spot.loc[spot.event == 'h', 'utc_time'].iloc[0] > sh.utc_time.iloc[0]:
            d = (spot.loc[spot.event == 'h', 'utc_time'].iloc[0] - sh.utc_time.iloc[0]).seconds
            y = 1 - 2*d / dh
        else:
            d = (sh.utc_time.iloc[0] - spot.loc[spot.event == 'h', 'utc_time'].iloc[0]).seconds
            y = -1 + 2*d / dh
        if not rising:
            y = -y

        x, y = rot2d((x, y), incl, deg=False)
        out.append({'name': name,
                    'x': x,
                    'y': y,
                    'incl_obs': None if complemented else (np.rad2deg(incl) if deg else incl),
                    'r_sun': r_sun,
                    'complemented': complemented})
    return out

def get_coords_rh(rec, alpha=None, r_sun=65, pole=None, deg=True):
    """Compute x, y coordinates from observations with rhomboid wires.
    The x-axis is aligned with the line of motion of the Sun, the y-axis points
    to the northern hemisphere, and the radius of the Sun is one.

    Parameters
    ----------
    rec : pandas.DataFrame
        Record with contact times.
    alpha : float
        Inclination of the side of the rhombus. Default arctan(2).
    r_sun : float
        Solar radius to be used for incomplete records.
    pole : str, optional
        Indicates the pole that contacts the horizontal wire. 'S' for the south pole,
        'N' for the north pole. If not given, the value from the record is used.
    deg : bool
        If True, angles are degrees, else radians. Default True.

    Returns
    -------
    out : list
        Dictionary with x, y coordinates and additional parameters for each sunspot in the record.
    """
    if alpha is None:
        alpha = np.arctan(2)
    elif deg:
        alpha = np.deg2rad(alpha)
    spot_names = list(rec.object.loc[rec.object != 'S'].unique())
    sun = rec.loc[rec.object == 'S']
    if pole is None:
        pole = rec.pole.iloc[0]
    sign = 1 if pole == 'S' else -1
    complemented = False
    err = 0

    times = sun.utc_time.values
    if not np.all(times[:-1] <= times[1:]):
        raise ValueError('Sun times order is broken in record {}.'.format(rec.obs_id.iloc[0]))

    if len(sun) == 4:
        d1 = (sun.utc_time.iloc[2] - sun.utc_time.iloc[0]).seconds
        d2 = (sun.utc_time.iloc[3] - sun.utc_time.iloc[1]).seconds
        err = np.arctan(np.tan(alpha) * (d1 - d2) / (d1 + d2))
        r_sun = d1 * np.sin(alpha - err) / 2
    elif len(sun) == 3:
        complemented = True
        d1 = (sun.utc_time.iloc[2] - sun.utc_time.iloc[1]).seconds
        d2 = d1
        r_sun = d1 * np.sin(alpha) / 2
        print('Warning: check sun times order in record {}.'.format(rec.obs_id.iloc[0]))
    elif len(sun) == 2:
        complemented = True
    else:
        raise ValueError('Incomplete record {}.'.format(rec.obs_id.iloc[0]))

    out = []
    for name in spot_names:
        spot = rec.loc[rec.object == name]
        if len(spot) != 2:
            print('Incomplete spot {} in record {}.'.format(name, rec.obs_id.iloc[0]))
            continue
        times = spot.utc_time.values
        if not np.all(times[:-1] <= times[1:]):
            raise ValueError('Spot {} times order is broken in record {}.'.format(name, rec.obs_id.iloc[0]))
        p1 = (spot.utc_time.iloc[0] - sun.utc_time.iloc[0]).seconds
        p2 = (sun.utc_time.iloc[-1] - spot.utc_time.iloc[-1]).seconds
        x = ((r_sun*(np.tan((alpha - err)/2)*np.tan(alpha - err) +
                     -np.tan((alpha + err)/2)*np.tan(alpha + err)) +
              p2*np.tan(alpha + err) - p1*np.tan(alpha - err)) /
             (np.tan(alpha - err) + np.tan(alpha + err)))
        y = sign * np.tan(alpha - err)*(x - r_sun*np.tan((alpha - err)/2) + p1)
        x /= r_sun
        y = y / r_sun - sign

        out.append({'name': name,
                    'x': x,
                    'y': y,
                    'r_sun': r_sun,
                    'err': np.rad2deg(err) if deg else err,
                    'complemented': complemented})
    return out

def get_solar_coords(rec, obs_loc, method, **kwargs):
    """Compute heliographic coordinates from observations with contact times.

    Parameters
    ----------
    rec : pandas.DataFrame
        Record with contact times.
    obs_loc : astropy.coordinates.earth.EarthLocation
        Observer's location.
    method : str
        Method of observations. Can be 'hv' for horizontal and vertical wires,
        'ob' for oblique wires, 'rh' for rhomboid wires.
    kwargs : dict, optional
        Additional keywords for the reconstuction method.

    Returns
    -------
    out : list
        Dictionary with x, y coordinates and additional parameters for each sunspot in the record.
    """
    if rec.empty:
        raise ValueError('Empty record.')
    if not len(rec.obs_id.unique()) == 1:
        raise ValueError('Expected a single record, found {}.'.format(len(rec.obs_id.unique())))
    i = rec.utc_time.argmin()
    dt = rec.utc_time.iloc[i]
    dt_rec = rec.recorded_time.iloc[i]
    b0, l0, p0 = B0(dt).deg, L0(dt).deg, P(dt).deg
    orient = orientation(obs_loc, dt).deg
    incl = orient + p0

    dt2 = rec.utc_time.max()
    incl2 = orientation(obs_loc, dt2).deg + P(dt2).deg
    incl_diff = abs(incl2 - incl)

    if method == 'ob':
        get_coords = get_coords_ob
    elif method == 'hv':
        get_coords = partial(get_coords_hv, incl=incl, deg=True)
    elif method == 'rh':
        get_coords = get_coords_rh
    else:
        raise ValueError('Unknown method {}.'.format(method))

    out = []
    for coords in get_coords(rec, **kwargs):
        vec = np.array((coords['x'], coords['y']))
        coords.update(xy_to_hc(vec, b0, l0, p0, deg=True))
        coords['l0'] = l0
        coords['b0'] = b0
        coords['p'] = p0
        coords['incl'] = incl
        coords['incl_diff'] = incl_diff
        coords['obs_id'] = rec.obs_id.iloc[0]
        coords['utc_time'] = dt
        coords['recorded_time'] = dt_rec
        coords['sidereal'] = rec.sidereal.iloc[0]
        coords['rising'] = rec.rising.iloc[0]
        coords['type'] = method
        coords['modified'] = ~rec.modified_time.isna().all()
        out.append(coords)
    return out

def make_grid(ax, step=15, c='gray', lw=0.5):
    """Make a grid on the solar disk."""
    s = 180 // step + 1
    phi = np.linspace(-np.pi/2, np.pi/2, 100)
    for t in np.linspace(-np.pi/2, np.pi/2, s):
        ax.plot(np.cos(t)*np.sin(phi), [np.sin(t)]*len(phi), c=c, lw=lw)
    t = np.linspace(-np.pi/2, np.pi/2, 100)
    for phi in np.linspace(-np.pi/2, np.pi/2, s):
        ax.plot(np.cos(t)*np.sin(phi), np.sin(t), c=c, lw=lw)

def make_disk(figsize=(7, 7), step=15, c='gray', lw=0.5):
    """Make solar disk."""
    _, ax = plt.subplots(figsize=figsize)
    disk = plt.Circle((0, 0), 1, color='#f5ed76')
    ax.add_patch(disk)
    if step is not None:
        make_grid(ax, step=step, c=c, lw=lw)
    ax.set_xticks([])
    ax.set_yticks([])
    return ax

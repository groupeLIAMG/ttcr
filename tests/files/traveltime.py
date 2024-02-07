# -*- coding: utf-8 -*-
"""
Traveltime:
Computing traveltimes in a homogeneous medium for a 3-term phase velocity.

@author: Björn Rommel
@email: rommel@seisrock.com
@version: 1.1.0
@date: 22.11.2023
"""


# --- import --- import --- import --- import --- import --- import --- import


import numpy as np                     # numpy library for numerics
from matplotlib import pyplot as plt   # pyplot graphics
from scipy import interpolate as si    # scipy.interpolate


# --- user-defined constants --- user-defined constants --- user-defined


# Define "Dog Creek Shale" (Thomsen, 1986)
TITLE = "Dog Creek Shale"   # title
VP0 = 1875.                 # reference P-velocity
VS0 = 826.                  # reference S-Velocity
DELTA = 0.100               # Thomsen's delta
EPSILON = 0.225             # Thomsen's epsilon


# Switch plots on/off
VELPLOT = True         # plot energy velocity
TIMEPLOT = True        # plot offset-time if more than 1 sample
SEG = True             # follow SEG guideline
FIGWIDTH = 3.537       # width of figure: 1 column of a Geophysics paper
FIGHEIGHT = FIGWIDTH   # width of figure


# Define the geometry source-to-receiver in a homogenous medium
# ... to create traveltime-offset plot
XXX1 = [0. + iii * 5 for iii in range(401)]   # equally spaced receivers
ZZZ1 = [1875.] * 401                          # at constant vertical distance
# ... for a single event
XXX2 = 1875.
ZZZ2 = 1875.


# Usage
# There is no easy way to compute the phase velocity for a given propagation
# direction of the energy velocity. Take a look at Berryman's (1979) paper: the
# phase angle would be the inverse of the tangent of some function of the phase
# velocity and its derivative! Instead, I pre-compute the phase velocity for
# some small increments in incidence angles, pre-compute the resultant energy
# velocity both in magnitude and propagation direction, and construct an
# interpolating spline polynomial. With help of that interpolating spline
# polynomial I can calculate the energy velocity for any given propagation
# direction. So, precompute the energy velocity from angle.START to angle.END
# and construct the scipy CubicSpline interpolation polynomial. Input are
# either reference velocities and Thomsen's anisotropy parameters,

# ###   scs = precompute(
# ###       vp0=VP0, vs0=VS0, delta=DELTA, epsilon=EPSILON, wavetype='P'),

# or the generic reference velocity and anisotropy parameters

# ###   scs = precompute(
# ###       v00=v00, rrr=rrr)

# where

# ###   v00, rrr = thomsen2generic(
# ###       vp0=VP0, vs0=VS0, delta=DELTA, epsilon=EPSILON, wavetype='P'),

# From there compute the traveltime,

# ###   time = traveltime(
# ###       xxx=XXX, zzz=ZZZ, scs=scs)

# You can replace wavetype='P' with wavetype='SV' for an SV-wave, but, here,
# Dog Creek Shale shows a false triplication that cannot yet be handled.


# Definition of Generic Anisotropy Parameters
# This program was written for use with ttcpry, and in ray-tracing we use the
# linearized phase velocity v = v0 (1 + r2 * sin(theta)^2 + r4 * sin(theta)*4).
# Note, this version differs from the one used in the paper; that one is based
# on the squared phase velocity. Hence. (r2, r4), here, are
# 1/2 * (r2, r4) of the paper!
# (https://github.com/groupeLIAMG/ttcr)


# --- do not change unless you know --- do not change unless you know --- do


# Incidence angle (for pre-computing energy velocities)
START = -180.   # first incidence angle (float)
NOS = 3601      # increment of incidence angles (int)
END = +180.     # last incidence angle (float)


# --- functions --- functions --- functions --- functions --- functions ---


def precompute(
        vp0=VP0, vs0=VS0, delta=DELTA, epsilon=EPSILON, wavetype='P',
        v00=VP0, rrr=None):
    """
    Compute an interpolation polynomial for the energy velocity.

    Inputs
    ------
        vp0, vs0, delta, epsilon, wavetype   # Thomsen's anisotropy parameters
        v00, rrr                             # generic anisotropy parameters

    Parameters
    ----------
    vp0 : float
        P-reference phase velocity. The default is VP0.
    vs0 : float
        S-reference phase velocity. The default is VS0.
    delta : float
        Thomsen's anisotropy parameter delta. The default is DELTA.
    epsilon : float
        Thomsen's anisotropy parameter epsilon. The default is EPSILON.
    v00 : float
        Generic reference velocity
    rrr : np.array ([5])
        Generic anisotropy parameters
    wavetype : char ('P' or 'SV')
        P- or SV wavetype requested. The default is a P-wave

    Returns
    -------
    scs : PPoly
        Cubic spline interpolation.
        See https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html

    """
    # pylint:disable=too-many-arguments
    # compute the series of incidence angles
    rad = angleseries(start=START, end=END, nos=NOS)
    # compute generic anisotropy parameters for a P-wave
    if rrr is None:
        v00, rrr = (
            thomsen2generic(
                vp0=vp0, vs0=vs0, delta=delta, epsilon=epsilon,
                wavetype=wavetype))
    # compute phase velocity
    phase = phasevelocity(v00=v00, rrr=rrr, rad=rad)
    # compute energy velocity
    energy = energyvelocity(phase=phase)
    # compute the interpolation polynom
    scs = si.CubicSpline(energy['rad'], energy['mag'], bc_type='periodic')
    # plot energy velocity
    if VELPLOT:
        polarplot(vel=energy)   # comment out if not desired
    # return interpolation polynom
    return scs


def angleseries(start=START, end=END, nos=NOS):
    """
    Compute a series of incidence angles.

    Parameters
    ----------
    start : float
        First incidence angle. The default is START.
    end : float
        Last incidence angle. The default is END.
    nos : int
        Number of incidence angles. The default is NOS.

    Returns
    -------
    rad : np.array ([nos-many])
        Array of incidence angles in rad.

    """
    # generate a series of incidence angles
    step = (end - start) / (nos - 1)                       # step size
    grad = (                                               # series in grad
        np.array(
            [start + step * iii for iii in range(nos)]))
    rad = np.deg2rad(grad)                                 # convert to radiant
    # return
    return rad


def thomsen2generic(
        vp0=VP0, vs0=VS0, delta=DELTA, epsilon=EPSILON, wavetype='P'):
    """
    Convert Thomsen's (1986) anisotropy parameters into generic ones.

    The definition of anisotropy parameters follows Thomsen, L., 1986, Weak
    elastic anisotropy: GEOPHYSICS, 51, 1954-1966.

    Note, a generic phase velocity shall be defined as
    vp(theta) = v0 (1 + r2 * sin(theta)^2 + r4 * sin(theta)^4)

    The generic anisotropy coefficients returned are rrr[0] = 1, rrr[2]=r2,
    rrr[4] = r4 to conform with Python notation.

    Parameters
    ----------
    vp0 : float
        P-reference phase velocity. The default is VP0.
    vs0 : float
        S-reference phase velocity. The default is VS0.
    delta : float
        Thomsen's anisotropy parameter delta. The default is DELTA.
    epsilon : float
        Thomsen's anisotropy parameter epsilon. The default is EPSILON.
    wavetype : char ('P' or 'SV')
        P- or SV wavetype requested. The default is a P-wave

    Returns
    -------
    v00 : float
        reference velocity
    rrr : np.array ([5])
        generic anisotropy coefficients for either a P- or SV-wave

    """
    # init
    rrr = np.full(5, fill_value=np.NAN, dtype=float)
    rrr[0] = 1
    # shortcut
    secterm = 1 + 2 * vp0 ** 2 / (vp0 ** 2 - vs0 ** 2) * delta
    # for a P-wave
    if wavetype == 'P':
        v00 = vp0
        rrr[2] = delta
        rrr[4] = (epsilon - delta) * secterm
    # for an SV-wave
    if wavetype == 'SV':
        v00 = vs0
        rrr[2] = vp0 ** 2 / vs0 ** 2 * (epsilon - delta)
        rrr[4] = -1 * rrr[2] * secterm
    # return
    return v00, rrr


def phasevelocity(v00=VP0, rrr=None, rad=None):
    """
    Compute the magnitude and incidence angle of the phase velocity.

    Parameters
    ----------
    v00 : float
        Reference velocity. The default is VP0.
    rrr : np.array ([5])
        Generic anisotropy parameters. The default is None.
        rrr must be first computed with help of thomsen2generic().
    rad : np.array ([nos-many])
        Incidence angles in rad. The default is None.
        rad must first have been computed with help of angle().

    Returns
    -------
    phase : dict
        Phase velocity with
        'mag' : absolute magnitude of the velocity
        'dmag' : first derivative of the magnitude with respect to the
            incidence angles (for use with the energy velocity)
        'rad' : incidence angle(s)
        'xxx' : horizontal component of the phase velocity
        'zzz' : vertical component of the phase velocity

    """
    # shortcuts
    sin = np.sin(rad)
    cos = np.cos(rad)
    # compute phase velocity
    phase = {
        # magnitude
        'mag': v00 * (1. + rrr[2] * sin ** 2 + rrr[4] * sin ** 4),
        # first derivative of the magnitude with respect to the incidence angle
        'dmag' : 2 * v00 * sin * cos * (rrr[2] + 2 * rrr[4] * sin ** 2),
        # direction of propagation
        'rad': rad}
    # compute components
    phase.update({
        'xxx': phase['mag'] * sin,
        'zzz': phase['mag'] * cos})
    # return
    return phase


def energyvelocity(phase=None):
    """
    Compute the magnitude and incidence angle of the energy velocity.

    The calculation follows Berryman, J. G., 1979, Long-wave elastic anisotropy
    in transversely isotropic media: GEOPHYSICS, 44, 896-917.

    Parameters
    ----------
    phase : dict
        Phase velocity. The default is None.
        phase must first be computed with help of phasevelocity().

    Returns
    -------
    energy : dict
        Energy velocity with
        'mag' : absolute magnitude of the velocity
        'rad' : direction of propagation or incidence angle(s)
        'xxx' : horizontal component of the energy velocity
        'zzz' : vertical component of the energy velocity

    """
    # compute magnitude and incidence angle of theenergy velocity
    tan = np.tan(phase['rad'])
    term = phase['dmag'] / phase['mag']
    energy = {
        # magnitude
        'mag': np.sqrt(phase['mag'] ** 2 + phase['dmag'] ** 2),
        # direction of propagation
        'rad': np.arctan((tan + term) / (1. - tan * term))
        }
    # correct energy angle
    energy['rad'] -= [
        np.pi if rad < -1. * np.pi / 2.
        else 0.
        for rad in phase['rad']]
    energy['rad'] += [
        np.pi if rad > +1. * np.pi / 2.
        else 0.
        for rad in phase['rad']]
    # compute components
    energy.update({
        'xxx': energy['mag'] * np.sin(energy['rad']),
        'zzz': energy['mag'] * np.cos(energy['rad'])})
    # return
    return energy


def polarplot(vel=None, title=TITLE):
    """
    Plot a polar diagram of the velocity.

    Parameters
    ----------
    vel : dict
        Phase or energy velocity. The default is None.
        Note, vel['zzz'] over vel['xxx'] will be plotted.
    title : char
        Title. The default is TITLE.

    Returns
    -------
    None.

    """
    # set up the figure
    fig = plt.figure(0, figsize=(FIGWIDTH, FIGWIDTH),layout='tight')
    axx = fig.add_subplot(111)
    # plot the velocity
    plt.plot(vel['xxx'], vel['zzz'])
    # find some reasonable velocity larger than the max actual one
    maxi = np.amax(np.sqrt(vel['xxx'] ** 2 + vel['zzz'] ** 2))
    maxi = np.ceil(maxi / 500) * 500
    # set axes limits
    plt.xlim(left=-1*maxi, right=+1*maxi)
    plt.ylim(bottom=-1.*maxi, top=+1*maxi)
    # set aspect ratio
    axx.set_aspect('equal', 'box')
    # set xlabel
    plt.xlabel("Horizontal Velocity")
    # set ylabel
    plt.ylabel("Vertical Velocity")
    # plot title
    plt.title(title + ": Polar Velocity")
    # save
    plt.savefig(
        'velocity.eps',
        format='eps', bbox_inches='tight', pad_inches=0., dpi=1200.)
    # show
    plt.show()

def traveltime(xxx=None, zzz=None, scs=None, seg=False):
    """
    Compute the traveltime for given horizontal / vertical distances.

    Parameters
    ----------
    xxx : list of float
        Horizontal distance from source to receiver(s). The default is None.
    zzz : list of float
        Vertical distance from source to receiver(s). The default is None.
    scs : PPoly
        Scipy cubic spline interpolation. The default is None.
    seg : boolean
        True / False for conforming with SEG guidelines.

    Returns
    -------
    time : (np.array of) float
        Traveltimes.

    """
    # convert to array
    xxx = np.array(xxx)
    zzz = np.array(zzz)
    # compute length and direction of the travel path
    length, angle = travelpath(xxx=xxx, zzz=zzz)
    # interpolate an energy velocity
    velocity = scs(angle, 0)
    # compute traveltime
    time = length / velocity
    # plot
    if TIMEPLOT:                               # check plotting switch
        if xxx.size > 1:                       # check sample number (> 1)
            travelplot(
                xxx=xxx, time=time, seg=seg)
    # return
    return time


def travelpath(xxx=None, zzz=None):
    """
    Compute the travel distance and direction of propagation.

    Parameters
    ----------
    xxx : np.array or list
        Horizontal distance from source to receiver. The default is None.
    zzz : np.array or list
        Vertical distance from source to receiver. The default is None.

    Returns
    -------
    length : np.array of float
        Length of travel path
    angle : np.array of float
        Direction of propagation (for which an energy velocity is requested)

    """
    # compute the length
    length = np.sqrt(xxx ** 2 + zzz ** 2)
    # compute the direction of propagation
    angle = np.arctan(xxx/zzz)
    # return
    return length, angle


def travelplot(xxx=None, time=None, title=TITLE, seg=False):
    """
    Plot traveltime over offset.

    Parameters
    ----------
    xxx : np.array of float
        Horizontal distance from source to receiver. The default is None.
    time : np.array of float
        Traveltime. The default is None.
    title : char
        Title. The default is TITLE.
    seg : boolean
        True / False for conforming with SEG guidelines.

    Returns
    -------
    None.

    """
    # pylint:disable=unused-variable
    # set up the figure
    fig = plt.figure(1, figsize=(FIGWIDTH, FIGHEIGHT), layout='tight')
    axx = fig.add_subplot(111)   # noqa   # let's have it defined
    # plot traveltime-offset curve
    plt.plot(xxx, time)
    # set offset limits
    maxi = np.ceil(xxx[-1] / 100) * 100
    plt.xlim(left=0., right=maxi)
    # set time limits
    maxi = np.ceil(time[-1] / 0.1) * 0.1
    plt.ylim(bottom=0., top=maxi)
    # place time axis
    plt.tick_params(
        axis='y', which='both', labelleft=True, labelright=True, right=True)
    # set offset label
    plt.xlabel("Offset")
    # set time label
    plt.ylabel("Traveltime")
    # plot title
    plt.title(title + ": Traveltime versus Offset")
    # SEG style
    if seg:
        setfont()
    # show
    plt.show()


def setfont():
    """
    Define the font styles, sizes to conform with SEG's guidelines.

    Returns
    -------
    None.

    """
    # set fonts https://stackoverflow.com/a/39566040
    plt.rc('text', usetex=True)          # Latex for versatility
    plt.rc('font', family='Helvetica')   # fontstyle
    plt.rc('ps', usedistiller='xpdf')    # avoiding bitmap
    plt.rc('font', size=8)               # controls default text sizes
    plt.rc('axes', titlesize=8)          # fontsize of the axes title
    plt.rc('axes', labelsize=8)          # fontsize of the axes labels
    plt.rc('xtick', labelsize=8)         # fontsize of the tick labels
    plt.rc('ytick', labelsize=8)         # fontsize of the tick labels
    plt.rc('legend', fontsize=8)         # legend fontsize
    plt.rc('figure', titlesize=8)        # fontsize of the figure title


# --- main --- main --- main --- main --- main --- main --- main --- main ---


def example1():
    """
    Demo the use of this module for use with Thomsen's anisotropy parameters.

    """
    # precompute the energy velocity and
    # construct the scipy CubicSpline interpolation polynomial
    # 1. starting from Thomsen's anisotropy parameters
    scs = (                                   # precompute the interpolation
        precompute(                           # polynomial from Thomsen's
            vp0=VP0, vs0=VS0, delta=DELTA,    # anisotropy parameters
            epsilon=EPSILON, wavetype='P'))
    # compute traveltime(s)
    time = (
        traveltime(
            xxx=XXX1, zzz=ZZZ1, scs=scs, seg=SEG))
    # return
    return time


def example2():
    """
    Demo the use of this module for use while splitting the precomputation.

    Use this example to extract generic anisotropy parameters.

    """
    # precompute the energy velocity and
    # construct the scipy CubicSpline interpolation polynomial
    # 2. starting from Thomsen's anisotropy parameters, but extracting generic
    #    anisotropy parameters in between
    v00, rrr = (                              # compute generic anisotropy
        thomsen2generic(                      # parameters from Thomsen's one
            vp0=VP0, vs0=VS0, delta=DELTA,    # or skip this step if known
            epsilon=EPSILON, wavetype='P'))
    print(f"\ndelta = {DELTA},   epsilon = {EPSILON}")
    print(f"\nrrr = {rrr}")
    scs = (                      # continue to precompute the
        precompute(              # interpolation polynomial from
            v00=v00, rrr=rrr))   # generic anisotropy parameters
    # compute traveltime(s)
    time = (
        traveltime(
            xxx=XXX2, zzz=ZZZ2, scs=scs, seg=SEG))
    # return
    return time


def example3():
    """
    Demo the use of this module for use with generic anisotropy parameters.

    Below, omit the call to thomsen2generic if the generic anisotropy
    parameters are known. Recall, all information about the wavetype is lost
    after that.

    """
    # precompute the energy velocity and
    # construct the scipy CubicSpline interpolation polynomial
    # 3. starting directly from generic anisotropy parameters
    v00 = VP0
    rrr = np.array([1., np.NAN, 0.1, np.NAN, 0.15602005])
    scs = (                      # continue to precompute the
        precompute(              # interpolation polynomial from
            v00=v00, rrr=rrr))   # generic anisotropy parameters
    # compute traveltime(s)
    time = (
        traveltime(
            xxx=XXX1, zzz=ZZZ1, scs=scs, seg=SEG))
    # return
    return time


if __name__ == "__main__":
    # run example 1
    print("\nrunning example 1:")
    time1 = example1()
    print(f"traveltime at {XXX1[375]}m offset: {time1[375]}")
    # run example 2
    print("\nrunning example 2:")
    time2 = example2()
    print(f"traveltime at {XXX2}m offset: {time2}")
    # run example 3 (most simila to ttcrpy)
    print("\nrunning example 3:")
    time3 = example3()
    print(f"traveltime at {XXX1[375]}m offset: {time3[375]}")
    print(f"traveltime at {XXX1[400]}m offset: {time3[400]}")

# -*- coding: utf-8 -*-
"""
Computing linearized traveltimes in a homogeneous medium.

This program was written for use with ttcpry, software for ray-tracing
simulation https://ttcrpy.readthedocs.io/en/latest/ of wavefield propagation,
by Berard Giroux https://inrs.ca/en/research/professors/bernard-giroux/

Here, we use the linearized phase velocity v = v0 (1 + r2 * sin(theta)^2 + r4 *
sin(theta)*4), where r2, r4 are generic anisotropy parameters. Note, this
version differs from the one used in the paper: (r2, r4), here, are
1/2 * (r2, r4) of the paper!

The software repository is https://github.com/groupeLIAMG/ttcr


@author: Björn Rommel
@email: seisrock@seisrock.com
@version: 1.0.0
@date: 28.11.2023
"""


# --- import --- import --- import --- import --- import --- import --- import


import numpy as np                 # numpy library for numerics
from scipy import optimize as so   # optimization library
from scipy import interpolate as si   # interpolation library


# --- functions --- functions --- functions --- functions nctions --- functions ---


def thomsen2generic(
        vp0=None, vs0=None, delta=None, epsilon=None, wavetype=None):
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


def energyvelocity(v00=None, rrr=None, theta=None):
    """
    Compute the magnitude and incidence angle of the phase velocity.

    Parameters
    ----------
    v00 : float
        Reference velocity. The default is VP0.
    rrr : np.array ([5])
        Generic anisotropy parameters. The default is None.
        rrr must be first computed with help of thomsen2generic().
    theta : np.array ([nos-many])
        Incidence angles in rad. The default is None.

    Returns
    -------
    phase : dict
        Phase velocity with
        'mag' : absolute magnitude of the velocity
        'dmag' : first derivative of the magnitude with respect to the
            incidence angles (for use with the energy velocity)
        'theta' : incidence angle(s)
        'xxx' : horizontal component of the phase velocity
        'zzz' : vertical component of the phase velocity

    """
    # shortcuts
    sin = np.sin(theta)
    # compute energy velocity
    velocity = v00 * (1. + rrr[2] * sin ** 2 + rrr[4] * sin ** 4)
    # return
    return velocity


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
    theta : np.array of float
        Direction of propagation (for which an energy velocity is requested)

    """
    # convert to np.array
    xxx = np.array(xxx)
    zzz = np.array(zzz)
    # compute the length
    length = np.sqrt(xxx ** 2 + zzz ** 2)
    # compute the direction of propagation
#    theta = np.arctan(xxx/zzz)
    theta = np.arctan2(xxx, zzz)
    # return
    return length, theta


def root(var, phi, rrr):
    """
    Compute the objective function for the root-finding algorithm.

    Parameters
    ----------
    var : np.array
        Current array of probing phase angles.
    phi : np.array
        Given array of energy angles.
    rrr : list
        Anisotropy parameters, see thomsen2generic.

    Returns
    -------
    ret : np.array
        Difference between given and probing energy angles.

    """
    # probing energy angle
    val = (
        np.arctan(
            np.tan(var) * (1 + 2 * rrr[2] + 4 * rrr[4] * np.sin(var) ** 2)))
    # correct into the same quadrant as the given energy angles
    iii = phi < -1 * np.pi / 2.
    val[iii] -= np.pi
    iii = phi > +1 * np.pi / 2.
    val[iii] += np.pi
    # difference
    ret = phi - val
    # return difference
    return ret


def mysolve(phi, rrr):
    """
    Invert the energy angle for the phase angle via interpolation.

    Parameters
    ----------
    phi : np.array
        Energy angles
    rrr : list
        Anisotropy parameters, see thomsen2generic.

    Returns
    -------
    theta : np.array
        Traveltime(s).

    """
    # set phase angle
    imax = 181   # number of samples (arbitrary positive integer)
    var = (      # increment from -pi/2 to ++pi/2 in steps of pi/(imax-1))
        np.array([
            -1 * np.pi / 2. + np.pi * iii / (imax - 1)
            for iii in range(imax)]))
    # compute energy angle
    val = (
        np.arctan(
            np.tan(var) * (1 + 2 * rrr[2] + 4 * rrr[4] * np.sin(var) ** 2)))
    # interpolate
    scs = si.CubicSpline(val, var)
    theta = scs(phi, 0)
    # return
    return theta


def ttwean(
        vp0=None, vs0=None, delta=None, epsilon=None, wavetype=None, xxx=None,
        zzz=None):
    """
    Invert an energy travelpath into a phase angle.

    The first algorithm used is basically a root-finding one: minimize the 
    difference between the given energy angle and the energy angle computed 
    from a probing phase angle.
    The second algorithm computes a semi-circle of energy angles versus phase
    angles, constructs an spline, and interpolates the phase angles for the
    given energy angles. This algorithm can be outsourced and pre-computed.

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
        P- or SV wavetype requested. The default is a P-wave.

    """
    # pylint:disable=too-many-arguments
    # compute reference velocity and generic anisotropy parameters
    v00, rrr = (
        thomsen2generic(
            vp0=vp0, vs0=vs0, delta=delta, epsilon=epsilon, wavetype=wavetype))
    length, phi = travelpath(xxx=xxx, zzz=zzz)
    # compute phase angles corresponding to raypath segments
    theta = so.fsolve(root, phi, (phi, rrr))   # option 1, comment out 1 or 2
    theta = mysolve(phi, rrr)                  # option 2, comment out 1 or 2
    # compute linearized energy velocities
    velocity = energyvelocity(v00=v00, rrr=rrr, theta=theta)
    # compute traveltimes
    time = length / velocity
    # return
    return time


# -----------------------------------------------------------------------------


def example1():
    """
    Demo the computation of a traveltime over a segment of a ray path.

    Returns
    -------
    None.

    """
    # Define "Dog Creek Shale" (Thomsen, 1986)
    vp0 = 1875.       # reference P-velocity
    vs0 = 826.        # reference S-Velocity
    delta = 0.100     # Thomsen's delta
    epsilon = 0.225   # Thomsen's epsilon
    # Set coordinates
    xxx = 400.   # horizontal distance of ray segment (float, list or array)
    zzz = 500.   # vertical distance of ray segment (float, list or array)
    # Define wavetype
    wavetype = 'P'
    # You can replace wavetype='P' with wavetype='SV' for an SV-wave, but Dog
    # Creek Shale shows a false triplication at large offset that cannot be
    # handled by a linearized energy velocity.
    # Compute traveltime
    time = (
        ttwean(
            vp0=vp0, vs0=vs0, delta=delta, epsilon=epsilon, wavetype=wavetype,
            xxx=xxx, zzz=zzz))
    # return
    return time

def example2():
    """
    Demo the computation of a traveltime over a segment of a ray path.

    Returns
    -------
    None.

    """
    # Stretch "Dog Creek Shale" from above
    # The stretch factor is g=0.201216, which makes the stretched medium
    # almost "nearly-isotropic." The reason being, the virtual source (base
    # point of a ray segment) moves from gridpoint 500m to exactly 548m, and
    # a grid spacing of 1m marks the limit of what my hardware can handle.
    vp0 = 2055.          # reference P-velocity
    vs0 = 905.296        # reference S-Velocity
    delta = -0.000506    # Thomsen's delta
    epsilon = 0.107758   # Thomsen's epsilon
    # Set coordinates
    xxx = 400.   # horizontal distance of ray segment (float, list or array)
    zzz = 548.   # vertical distance of ray segment (float, list or array)
    # Define wavetype
    wavetype = 'P'
    # You can replace wavetype='P' with wavetype='SV' for an SV-wave, but Dog
    # Creek Shale shows a false triplication at large offset that cannot be
    # handled by a linearized energy velocity.
    # Compute traveltime
    time = (
        ttwean(
            vp0=vp0, vs0=vs0, delta=delta, epsilon=epsilon, wavetype=wavetype,
            xxx=xxx, zzz=zzz))
    # return
    return time


# -----------------------------------------------------------------------------


if __name__ == "__main__":
    # call example
    time1 = example1()
    time2 = example2()
    print(f"original traveltime:  {time1}")
    print(f"stretched traveltime: {time2}")

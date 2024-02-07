# -*- coding: utf-8 -*-
"""
ttwean: Velocity-versus-Structure for Reflection Seismics!

ttwean is a collection of function used to investigate "stretch", a method to
generate a series of velocity-structure models all of which (approximately)
satisfy the same traveltime curve.

This program was written for use with ttcpry, a software for ray-tracing
simulation https://ttcrpy.readthedocs.io/en/latest/ of wavefield propagation
by Berard Giroux https://inrs.ca/en/research/professors/bernard-giroux/
The software repository is https://github.com/groupeLIAMG/ttcr

You can use the functions ttinit and ttvel as part of the entire ttwean
package, which can be installed via pip ttwean...whl. However, for use
specifically with ttcrpy only you can safely rip out the functions ttinit and
ttvel below and use them independently outside that package.

Note, due to some misunderstanding of ttcrpy on my side, there is a conflict of
notation: the parameter r in ttcpry, see the example 5,
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example5.ipynb
ought to be the linear energy-velocity parameter s of my paper; the latter is
computed below and must be used in place of the parameter r in ttcpry!

My software repository is https://github.com/bjornrommel/steinkauz under the
project Subsurface-Velocity Ambiguity.

@author: Björn E. Rommel
@email: ttwean@seisrock.com
@version: 3.0.1
@date: 2024-02-03
"""


# ttinit   initialize energy-velocity parameters for each medium, wavetype, or
#          stretch; it is independent of any incidence angle or time step
# ttvel    compute the energy velocity for each medium, wavetype, or incidence
#          angle; in practice, probably for each cell, incidence angle and time
#          step

# Python imports
import numpy as np   # numpy


# define "Dog Creek Shale" (Thomsen, 1986)
# (any example will do, but need some pre-defined input)
TITLE = "Dog Creek Shale"   # title
VP0 = 1875.                 # reference P-velocity
VS0 = 826.                  # reference S-Velocity
DELTA = 0.100               # Thomsen's delta
EPSILON = 0.225             # Thomsen's epsilon
GGG = 0.201216              # stretch parameter used in my paper
GGG = 2 * DELTA             # stretch parameter to near isotropy for "P"
# ### GGG = 0.                    # neutral stretch parameter
WAVETYPE = "P"              # "P"-wave
# ### WAVETYPE = "SV"             # "SV"-wave
ANGLE = 30.                 # incidence angle (phase=energy angle)


class TTWean():
    """
    Characterize energy velocity for simulation with ttcpry.

    """

    # initialize
    def __init__(self):
        pass

    # set the energy-velocity parameters
    @staticmethod
    def ttinit(
            vp0:float=None, vs0:float=None, delta:float=0.,
            epsilon:float=0., wavetype:str="P", ggg:float=0.) -> np.array:
        """
        Precompute the linear energy-velocity parameters once.

        ttinit must be run for each medium, wavetype, and/or stretch once. It
        is the user's responsibility to assign the correct output to each cell.
        See Giroux's Jupyter notebook for details.
        ttvel must be run for each incidence angle, and it would be valid,
        albeit impractical, for all cells of the same medium and for the same
        wavetype and stretch. Typically, though, it would be run for each cell
        for each time step.

        Typically,
            vp0[medium,stretch], vs0[medium,stretch],
            sss[medium,stretch,wavetype] = (
                ttinit(
                    vp0=VP0, vs0=VS0, delta=DELTA, epsilon=EPSILON,
                    wavetype="P" or wavetype="SV", ggg=GGG))
            vel[wavetype,angle,cell] = (
                ttvel(
                    vp0=vp0[medium,stretch], vs0=vs0[medium,stretch],
                    sss=sss[medium,stretch], wavetype=wavetype, angle=ANGLE)
        Note, all capital for user-defined parameters.
        Any stretch modifies P and SV-velocities: so, be careful not to
        overwrite your in- and output variable.
        Use the same wavetype for ttinit and ttray.

        Parameters
        ----------
        vp0 : float
            reference P-velocity
        vs0 : float
            reference S-velocity
        delta : float, default is 0.
            Thomsen parameter delta
        epsilon : float, default is 0.
            Thomsen parameter epsilon
        wavetype : char
            "P" : P-wave
            "SV" : S-wave
        ggg : float, default is 0. (neutral)
            stretch factor

        Returns
        -------
        vp0, vs0 : float
            P- and SV-reference velocities: modified only if stretched
        sss : np.array ([5x1])
            linear energy-velocity parameters

        """
        # pylint:disable=too-many-arguments
        # shortcuts for computing squared phase-velocity parameters
        vps2 = (vp0 / vs0) ** 2
        fac = 1. + 2. * vps2 / (vps2 - 1.) * delta
        # define squared phase-velocity parameters
        rrr =  np.full((5,1), fill_value=np.NAN, dtype=float)
        rrr[0,0] = 1.
        # compute squared phase-velocity parameters for a qP-wave
        if wavetype == "P":
            rrr[2,0] = 2. * delta
            rrr[4,0] = 2. * (epsilon - delta) * fac
        # compute squared phase-velocity parameters for a qSV-wave
        if wavetype == "SV":
            rrr[2,0] = 2. * vps2 * (epsilon - delta)
            rrr[4,0] = -1. * rrr[2,0] * fac
        # stretch
        if ggg != 0.:
            # stretch reference velocity
            vp0 *= np.sqrt(1. + ggg)
            vs0 *= np.sqrt(1. + ggg)
            # stretch squared phase-velocity parameters
            rrr[2,0] = (rrr[2,0] - ggg) / (1. + ggg)
            rrr[4,0] = rrr[4,0] / ((1. + ggg) ** 2)
        # compute squared energy-velocity parameters
        ttt =  np.full((5,1), fill_value=np.NAN, dtype=float)
        ttt[0,0] = 1.
        ttt[2,0] = rrr[2,0] / (1. + rrr[2,0])
        ttt[4,0] = (
            (rrr[2,0] ** 2 * (1. + rrr[2,0]) ** 2 + rrr[4,0])
            /
            (1. + rrr[2,0]) ** 4)
        # compute linear energy velocity parameters
        sss =  np.full((5,1), fill_value=np.NAN, dtype=float)
        sss[0,0] = 1.
        sss[2,0] = ttt[2,0] / 2.
        sss[4,0] = -1. * ttt[2,0] ** 2 / 8. + ttt[4,0] / 2.
        # return
        return vp0, vs0, sss

    # compute energy velocity
    @staticmethod
    def ttvel(
            vp0:float=None, vs0:float=None,
            sss:np.array=np.array([[1.],[np.nan],[0.],[np.nan],[0]]),
            wavetype:str="P", angle:float=None) -> float:
        """
        Compute a linear energy velocity for a given directional angle.

        ttinit must be run for each medium, wavetype, and/or stretch once. It
        is the user's responsibility to assign the correct output to each cell.
        See Giroux's Jupyter notebook for details.
        ttvel must be run for each incidence angle, and it would be valid,
        albeit impractical, for all cells of the same medium and for the same
        wavetype and stretch. Typically, though, it would be run for each cell
        for each time step.

        Typically,
            vp0[medium,stretch], vs0[medium,stretch],
            sss[medium,stretch,wavetype] = (
                ttinit(
                    vp0=VP0, vs0=VS0, delta=DELTA, epsilon=EPSILON,
                    wavetype="P" or wavetype="SV", ggg=GGG))
            vel[wavetype,angle,cell] = (
                ttvel(
                    vp0=vp0[medium,stretch], vs0=vs0[medium,stretch],
                    sss=sss[medium,stretch], wavetype=wavetype, angle=ANGLE)
        Note, all capital for user-defined parameters.
        Any stretch modifies P and SV-velocities: so, be careful not to
        overwrite your in- and output variable.
        Use the same wavetype for ttinit and ttray.

        Parameters
        ----------
        vp0 : float
            reference P-velocity
        vs0 : float
            reference S-velocity
        sss : np.array ([5x1])
            linear energy-velocity parameters, default to isotropy
            note, s0=1 set, s2 and s4 in use, s1 and s3 nan
        wavetype : char
            "P" : P-wave
            "SV" : S-wave
        angle : float
            incidence angle in radiant

        Returns
        -------
        vel : float
            magnitude of an energy velocity for the given incidence angle,
            wavetype, and linear energy-velocity parameters.

        """
        # trig values for phase angles
        sin = np.sin(angle)
        sin2 = sin * sin
        # velocity
        vel = (
            {"P": vp0, "SV": vs0}[wavetype]               # reference velocity
            *                                             # *
            (1. + (sss[2,0] + sss[4,0] * sin2) * sin2))   # (angle-dependent)
        # return
        return vel


def main():
    """
    Entry point for demonstrating ttinit and ttvel for use with ttcrpy.

    Returns
    -------
    None.

    """
    # forward to pre-compute the linear energy-velocity parameters for each
    # medium, wavetype or stretch
    vp0, vs0, sss = (
        TTWean.ttinit(
            vp0=VP0, vs0=VS0, delta=DELTA, epsilon=EPSILON, wavetype=WAVETYPE,
            ggg=GGG))
    # forward to compute the linear energy velocity for each pre-/post-stretch
    # medium, wavetype or incidence angle
    vel = (
        TTWean.ttvel(
            vp0=vp0, vs0=vs0,# sss=sss, wavetype=WAVETYPE,
            angle=np.deg2rad(ANGLE)))
    # print output
    text = f"{TITLE}\n"                        # title
    text += f"post-stretch with g = {GGG}\n"   # stretch factor
    text += f"vp0 = {vp0}\n"                   # P-reference velocity
    text += f"vs0 = {vs0}\n"                   # S-reference velocity
    text += f"s\u2082 = {sss[2,0]}\n"          # lin. energy-velocity para. 2
    text += f"s\u2084 = {sss[4,0]}\n"          # lin. energy-velocity para. 4
    text += f"v({ANGLE}\u00b0) = {vel}"        # energy velocity (angle)
    print(text)


# --- main --------------------------------------------------------------------


if __name__ == "__main__":
    main()

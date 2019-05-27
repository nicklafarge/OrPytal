########### Standard ###########
import logging
import itertools
import unittest

########### Local ###########
from orpytal import units, Orbit, KeplarianState, bodies, plotting, get_plot_utils, integration, frames, Trajectory
from orpytal.errors import InvalidInputError

########### External ###########
import numpy as np

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)

earth = bodies.earth
orbit = Orbit(earth,
              a=51000 * units.km,
              e=0.7,
              raan=10 * units.deg,
              arg_periapsis=10 * units.deg,
              inclination=45 * units.deg)

st = orbit.get_state(ta=0.2)

print(st.position.inertial(st))
print(st.position.perifocal(st).inertial(st))

print(st.position.perifocal(st))
print(st.position.inertial(st).perifocal(st))


r = st.position.value
C_ri = frames.RotatingFrame.inertial_dcm(orbit, st)
C_ir = C_ri.transpose()

C_rp = frames.RotatingFrame.perifocal_dcm(orbit, st)
C_pr = C_rp.transpose()


C_pi = np.dot(C_pr, C_ri)
C_pi = np.dot(C_ri, C_pr)

ri = np.dot(C_ri, r)
rp = np.dot(C_rp, r)

ri_test = np.dot(C_pi, rp)

print("--actual values--")
print(ri)
print(ri_test)

print(rp)


def test_not_ascending(self):
    """"""
    orbit = Orbit(earth, a=a, e=0.4)
    st = orbit.get_state(r=49000)
    assert st._ascending == None
    assert st.ta == None

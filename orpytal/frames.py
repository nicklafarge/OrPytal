########### Standard ###########

########### External ###########
import numpy as np

from orpytal.errors import ParameterUnavailableError

########### Local ###########
from orpytal.utils import conics_utils


class Vector(object):
    def __init__(self,value, frame):
        self.value = value
        self.frame = frame

    @classmethod
    def from_vector(cls, orbit, state, vector):
        return cls(orbit, state, vector.value, vector.frame)

    def __str__(self):
        return '{} in {}'.format(self.value, self.frame.name)

    def __getitem__(self, i):
        return self.value.__getitem__(i)

    def __len__(self):
        return len(self.value)

    def __iter__(self):
        return self.value.__iter__()

    def __eq__(self, other):
        if isinstance(other, Vector):

            vec1 = self
            vec2 = other

            if vec1.frame != vec2.frame:
                return False
            else:
                return np.allclose(vec1.value, vec2.value)

        else:
            return False

    def __neq__(self, other):
        if isinstance(other, Vector):
            return not self == other
        return NotImplemented

    def to(self, frame, state):
        if frame == InertialFrame:
            dcm = self.frame.inertial_dcm
        elif frame == RotatingFrame:
            dcm = self.frame.rotating_dcm
        elif frame == PerifocalFrame:
            dcm = self.frame.perifocal_dcm
        else:
            raise ValueError('Frame ({}) is not recognized'.format(frame))

        if self.frame == frame:
            value = self.value
        else:
            value = dcm(state).dot(self.value.m) * self.value.units

        return Vector(value, frame)

    def inertial(self, *args):
        return self.to(InertialFrame, *args)

    def rotating(self, *args):
        return self.to(RotatingFrame, *args)

    def perifocal(self, *args):
        return self.to(PerifocalFrame, *args)

    def norm(self):
        return np.linalg.norm(self.value.m) * self.value.units

    def unit(self):
        return (self.value / np.linalg.norm(self.value.m)).m

    def cross(self, vector2):
        assert self.frame == vector2.frame
        units = self.value.units * vector2.value.units
        return np.cross(self.value, vector2.value)

    def dot(self, vector2):
        if isinstance(vector2, Vector):
            assert self.frame == vector2.frame
            units = self.value.units * vector2.value.units
            val2 = vector2.value
        else:
            units = 1
            val2 = vector2

        return np.dot(self.value.m, val2.m) * units


class CoordinateFrame(object):
    name = None
    fn_name = None

    @classmethod
    def inertial_dcm(cls, state):
        raise NotImplementedError()

    @classmethod
    def rotating_dcm(cls, state):
        raise NotImplementedError()

    @classmethod
    def perifocal_dcm(cls, state):
        raise NotImplementedError()

    @classmethod
    def vnc_dcm(cls, state):
        raise NotImplementedError()


class RotatingFrame(CoordinateFrame):
    name = 'Rotating Frame'
    fn_name = 'rotating'

    @classmethod
    def inertial_dcm(cls, state):
        requirements = ['raan', 'inclination', 'arg_latitude']
        if not conics_utils.state_orbit_satisfied(state, requirements):
            raise ParameterUnavailableError(
                'Need ascending node, inclination and argument of latitude to convert to xyz')

        Omega = state.orbit.raan
        theta = state.arg_latitude
        i = state.orbit.inclination
        dcm_ri = np.zeros((3, 3))
        dcm_ri[0, 0] = np.cos(Omega) * np.cos(theta) - np.sin(Omega) * np.cos(i) * np.sin(theta)
        dcm_ri[0, 1] = -np.cos(Omega) * np.sin(theta) - np.sin(Omega) * np.cos(i) * np.cos(theta)
        dcm_ri[0, 2] = np.sin(Omega) * np.sin(i)
        dcm_ri[1, 0] = np.sin(Omega) * np.cos(theta) + np.cos(Omega) * np.cos(i) * np.sin(theta)
        dcm_ri[1, 1] = -np.sin(Omega) * np.sin(theta) + np.cos(Omega) * np.cos(i) * np.cos(theta)
        dcm_ri[1, 2] = -np.cos(Omega) * np.sin(i)
        dcm_ri[2, 0] = np.sin(i) * np.sin(theta)
        dcm_ri[2, 1] = np.sin(i) * np.cos(theta)
        dcm_ri[2, 2] = np.cos(i)

        return dcm_ri

    @classmethod
    def rotating_dcm(cls, state):
        return np.eye(3)

    @classmethod
    def perifocal_dcm(cls, state):
        return PerifocalFrame.rotating_dcm(state).transpose()

    @classmethod
    def vnc_dcm(cls, state):
        return np.array(
            [[np.sin(state.fpa), 0, np.cos(state.fpa)],
             [np.cos(state.fpa), 0, -np.sin(state.fpa)],
             [0, 1, 0]])


class VncFrame(CoordinateFrame):
    name = 'Vnc Frame'
    fn_name = 'vnc'

    @classmethod
    def inertial_dcm(cls, state):
        dcm_vr = cls.rotating_dcm(state)
        dcm_ri = RotatingFrame.inertial_dcm(state)
        return dcm_vr.dot(dcm_ri)

    @classmethod
    def rotating_dcm(cls, state):
        return RotatingFrame.vnc_dcm(state).transpose()

    @classmethod
    def perifocal_dcm(cls, state):
        dcm_vr = cls.rotating_dcm(state)
        dcm_re = RotatingFrame.perifocal_dcm(state)
        return dcm_vr.dot(dcm_re)

    @classmethod
    def vnc_dcm(cls, state):
        return np.eye(3)


class PerifocalFrame(CoordinateFrame):
    name = 'Perifocal Frame'
    fn_name = 'perifocal'

    @classmethod
    def inertial_dcm(cls, state):
        dcm_pr = cls.rotating_dcm(state)
        dcm_ri = RotatingFrame.inertial_dcm(state)
        return dcm_ri.dot(dcm_pr) # order different from 440 because column vector

    @classmethod
    def rotating_dcm(cls, state):
        if state.ta is None:
            raise ParameterUnavailableError('Need true anomaly to convert to the rotating frame')
        return np.array(
            [[np.cos(state.ta).m, np.sin(state.ta).m, 0],
             [-np.sin(state.ta).m, np.cos(state.ta).m, 0],
             [0, 0, 1]]
        )

    @classmethod
    def perifocal_dcm(cls, state):
        return np.eye(3)

    @classmethod
    def vnc_dcm(cls, state):
        dcm_vr = cls.rotating_dcm(state)
        dcm_rv = RotatingFrame.vnc_dcm(state)
        return dcm_rv.dot(dcm_vr)


class InertialFrame(CoordinateFrame):
    name = 'Inertial Frame'
    fn_name = 'inertial'

    @classmethod
    def inertial_dcm(cls, state):
        return np.eye(3)

    @classmethod
    def rotating_dcm(cls, state):
        return RotatingFrame.inertial_dcm(state).transpose()

    @classmethod
    def perifocal_dcm(cls, state):
        return PerifocalFrame.inertial_dcm(state).transpose()

    @classmethod
    def vnc_dcm(cls, state):
        return VncFrame.inertial_dcm(state).transpose()

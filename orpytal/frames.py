########### Standard ###########

########### Local ###########
from orpytal.common import conics_utils
from orpytal.errors import ParameterUnavailableError

########### External ###########
import numpy as np


class Vector(object):
    def __init__(self, orbit, state, value, frame):
        self.orbit = orbit
        self.state = state

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
            while vec1.frame != vec2.frame:
                try:
                    vec1 = vec1.to(vec2.frame)
                    break
                except ParameterUnavailableError:
                    pass

                try:
                    vec2 = vec2.to(vec1.frame)
                    break
                except ParameterUnavailableError:
                    raise ParameterUnavailableError(
                        'Tried to Compare [{}] to [{}] Values, but could not put them in the same frame'.format(self, other))

            return np.allclose(vec1.value, vec2.value)

        return NotImplemented

    def __neq__(self, other):
        if isinstance(other, Vector):
            return not self == other
        return NotImplemented

    def to(self, frame):
        if frame == InertialFrame:
            dcm = self.frame.inertial_dcm
        elif frame == RotatingFrame:
            dcm = self.frame.rotating_dcm
        elif frame == OrbitFixedFrame:
            dcm = self.frame.orbit_fixed_dcm
        else:
            raise ValueError('Frame ({}) is not recognized'.format(frame))

        if self.frame == frame:
            value = self.value
        else:
            value = dcm(self.orbit, self.state).dot(self.value) * self.value.units

        return Vector(self.orbit, self.state, value, frame)

    def inertial(self):
        return self.to(InertialFrame)

    def rotating(self):
        return self.to(RotatingFrame)

    def orbit_fixed(self):
        return self.to(OrbitFixedFrame)

    def norm(self):
        return np.linalg.norm(self.value) * self.value.units

    def unit(self):
        return (self.value / np.linalg.norm(self.value)).m

    def cross(self, vector2):
        assert self.frame == vector2.frame
        units = self.value.units * vector2.value.units
        return np.cross(self.value, vector2.value) * units

    def dot(self, vector2):
        if isinstance(vector2, Vector):
            assert self.frame == vector2.frame
            units = self.value.units * vector2.value.units
            val2 = vector2.value
        else:
            units = 1
            val2 = vector2

        return np.dot(self.value, val2) * units



class CoordinateFrame(object):
    name = None
    fn_name = None

    @classmethod
    def inertial_dcm(cls, orbit, state):
        raise NotImplementedError()

    @classmethod
    def rotating_dcm(cls, orbit, state):
        raise NotImplementedError()

    @classmethod
    def orbit_fixed_dcm(cls, orbit, state):
        raise NotImplementedError()

    @classmethod
    def vnc_dcm(cls, orbit, state):
        raise NotImplementedError()

    def __mul__(self, other):
        return (self, other)

    # def __rmul__(self, other):
    #     self.__mul__(other)

class RotatingFrame(CoordinateFrame):
    name = 'Rotating Frame'
    fn_name = 'rotating'

    @classmethod
    def inertial_dcm(cls, orbit, state):

        requirements = ['ascending_node', 'inclination', 'arg_latitude']
        if not conics_utils.state_orbit_satisfied(state, orbit, requirements):
            raise ParameterUnavailableError('Need ascending node, inclination and argument of latitude to convert to xyz')

        Omega = orbit.ascending_node
        theta = state.arg_latitude
        i = orbit.inclination
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
    def rotating_dcm(cls, orbit, state):
        return np.eye(3)

    @classmethod
    def orbit_fixed_dcm(cls, orbit, state):
        return OrbitFixedFrame.rotating_dcm(orbit, state).transpose()

    @classmethod
    def vnc_dcm(cls, orbit, state):
        return np.array(
            [[np.sin(state.fpa), 0, np.cos(state.fpa)],
             [np.cos(state.fpa), 0, -np.sin(state.fpa)],
             [0, 1, 0]])


class VncFrame(CoordinateFrame):
    name = 'Vnc Frame'
    fn_name = 'vnc'

    @classmethod
    def inertial_dcm(cls, orbit, state):
        dcm_vr = cls.rotating_dcm(orbit, state)
        dcm_ri = RotatingFrame.inertial_dcm(orbit, state)
        return dcm_vr.dot(dcm_ri)

    @classmethod
    def rotating_dcm(cls, orbit, state):
        return RotatingFrame.vnc_dcm(orbit, state).transpose()

    @classmethod
    def orbit_fixed_dcm(cls, orbit, state):
        dcm_vr = cls.rotating_dcm(orbit, state)
        dcm_re = RotatingFrame.orbit_fixed_dcm(orbit, state)
        return dcm_vr.dot(dcm_re)

    @classmethod
    def vnc_dcm(cls, orbit, state):
        return np.eye(3)


class OrbitFixedFrame(CoordinateFrame):
    name = 'Orbit Fixed Frame'
    fn_name = 'orbit_fixed'

    @classmethod
    def inertial_dcm(cls, orbit, state):
        dcm_vr = cls.rotating_dcm(orbit, state)
        dcm_ri = RotatingFrame.inertial_dcm(orbit, state)
        return dcm_vr.dot(dcm_ri)

    @classmethod
    def rotating_dcm(cls, orbit, state):
        if state.ta is None:
            raise ParameterUnavailableError('Need true anomaly to convert to the rotating frame')
        return np.array(
            [[np.cos(state.ta), np.sin(state.ta), 0],
             [-np.sin(state.ta), np.cos(state.ta), 0],
             [0, 0, 1]]
        )

    @classmethod
    def orbit_fixed_dcm(cls, orbit, state):
        return np.eye(3)

    @classmethod
    def vnc_dcm(cls, orbit, state):
        dcm_vr = cls.rotating_dcm(orbit, state)
        dcm_rv = RotatingFrame.vnc_dcm(orbit, state)
        return dcm_vr.dot(dcm_rv)


class InertialFrame(CoordinateFrame):
    name = 'Inertial Frame'
    fn_name = 'inertial'

    @classmethod
    def inertial_dcm(cls, orbit, state):
        return np.eye(3)

    @classmethod
    def rotating_dcm(cls, orbit, state):
        return RotatingFrame.inertial_dcm(orbit, state).transpose()

    @classmethod
    def orbit_fixed_dcm(cls, orbit, state):
        return OrbitFixedFrame.inertial_dcm(orbit, state).transpose()

    @classmethod
    def vnc_dcm(cls, orbit, state):
        return VncFrame.inertial_dcm(orbit, state).transpose()

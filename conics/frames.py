########### Standard ###########

########### Local ###########

########### External ###########
import numpy as np


class Vector(object):
    def __init__(self, value, frame):
        self.value = value
        self.frame = frame

    def inertal(self):
        return self.frame.inertial_dcm.dot(self.value)

    def rotating(self):
        return self.frame.rotating_dcm.dot(self.value)

    def orbit_fixed(self):
        return self.frame.orbit_fixed_dcm.dot(self.value)


class CoordinateFrame(object):
    def __init__(self):
        pass

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


class RotatingFrame(CoordinateFrame):
    def __init__(self):
        super().__init__()

    @classmethod
    def inertial_dcm(cls, orbit, state):
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
    def __init__(self):
        super().__init__()

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
    def __init__(self):
        super().__init__()

    @classmethod
    def inertial_dcm(cls, orbit, state):
        dcm_vr = cls.rotating_dcm(orbit, state)
        dcm_ri = RotatingFrame.inertial_dcm(orbit, state)
        return dcm_vr.dot(dcm_ri)

    @classmethod
    def rotating_dcm(cls, orbit, state):
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
    def __init__(self):
        super().__init__()

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

########### Standard ###########
import math

########### Local ###########
from orpytal.common.utils import conics_utils

########### External ###########
import numpy as np
import scipy as sp


class KeplarianState(object):

    def __init__(self, orbit, name=''):
        self.orbit = orbit
        self.name = name

        self._ascending = None

        self._r = None
        self._v = None

        self._r_rth = None
        self._r_eph = None
        self._r_xyz = None

        self._v_rth = None
        self._v_eph = None
        self._v_xyz = None

        self._ta = None
        self._theta = None
        self._E = None
        self._M = None

        self._fpa = None
        self._ttp = None

        self._dcm_er = None
        self._dcm_rv = None
        self._dcm_ri = None

    def set_vars(self):
        self.r = None
        self.v = None

        self.r_rth = None
        self.r_eph = None
        self.r_xyz = None

        self.v_rth = None
        self.v_eph = None
        self.v_xyz = None

        self.ta = None
        self.theta = None
        self.E = None
        self.M = None
        self.fpa = None
        self.ttp = None

        self.dcm_er = None
        self.dcm_rv = None
        self.dcm_ri = None

        self.set_orbit()

    def set_orbit(self):
        self.orbit.p_from_state(self)
        self.orbit.e_from_state(self)
        self.orbit.se_from_state(self)
        self.orbit.h_from_state(self)
        self.orbit.h_xyz_from_state(self)

    @property
    def E(self):
        return self._E

    @E.setter
    def E(self, E=None):
        if self._E is not None or self.orbit.is_circular():
            return

        if E is not None:
            self._E = E
        elif self._M is not None:
            self._E = self._iteratively_find_E()
        elif self._ta is not None and self._r is not None and self.orbit.a is not None and self.orbit.e is not None:
            if self._ta == 0:
                self._E = 0
            elif self._ta == np.pi:
                self._E = np.pi
            else:
                self._E = self._get_ascending_sign* np.arccos(
                    (self.orbit.a - self._r) / (self.orbit.a * self.orbit.e))
        elif self.orbit.a is not None and self._r is not None and self.orbit.e is not None and self.orbit.is_ellipse():
            self._E = self._get_ascending_sign* np.arccos(((self.orbit.a - self._r) / (self.orbit.a * self.orbit.e)))
        else:
            return

        self.set_vars()

    def maneuver_vnc(self, delta_v):
        if self._dcm_rv is None:
            raise ValueError('Need DCM to Relate rth-vnc frames')

        print('Delta V (vnc): ' + str(delta_v))
        delta_v_rth = self._dcm_rv.dot(delta_v)
        return self.maneuver_rth(delta_v_rth)

    def maneuver_rth(self, delta_v):
        print('Delta V (rth): ' + str(delta_v))

        if self._dcm_er is None:
            raise ValueError('Need DCM to Relate eph-rth frames')

        delta_v_eph = self._dcm_er.transpose().dot(delta_v)

        return self.maneuver_eph(delta_v_eph)

    def maneuver_eph(self, delta_v):
        print('Delta V (eph): ' + str(delta_v))

        new_orbit = ConicOrbit(self.orbit.central_body)
        new_state = KeplarianState(new_orbit)
        new_state.r = self.r
        new_state.v = self.v - 0.5

        if self.dcm_ri is not None and self.dcm_er is not None:
            delta_v_xyz = self.dcm_ri.dot(self.dcm_er.dot(delta_v))
            new_state.v_xyz = self.v_xyz + delta_v_xyz
            new_state.r_xyz = self.r_xyz

            # TODO check sign here
            # new_state.set_fpa(self.fpa)
            # new_state.set_fpa(0)
            new_state.fpa = 0

            new_state.orbit.se_from_state(new_state)
            new_state.orbit.e_from_state(new_state)

            delta_omega = self._ta - new_state._ta
            new_state.orbit.arg_periapsis = self.orbit.arg_periapsis + delta_omega
            new_state.set_vars()
        return new_state

    def set_ascending(self):
        self._ascending = True

    def set_descending(self):
        self._ascending = False

    def _iteratively_find_E(self, tol=1e-12):
        return self._find_E_rec(self.M, self.M, tol=tol)

    def _find_E_rec(self, E, M, tol=1e-12):

        E1 = E - (E - self.orbit.e * math.sin(E) - M) / (1 - self.orbit.e * math.cos(E))
        dE = math.fabs(E - E1)

        if dE < tol:
            return E

        return self._find_E_rec(E1, M, tol=tol)

    def _is_ascending(self):
        if self._ascending is not None:
            return self._ascending

        angle_to_check = None

        if self._ta is not None:
            angle_to_check = self._ta

        elif self._E is not None:
            angle_to_check = self._E

        elif self._fpa is not None:
            angle_to_check = self._fpa

        return True if angle_to_check and 0 <= angle_to_check <= np.pi else False

    def _get_ascending_sign(self):
        return 1 if self._is_ascending() else -1

    def __repr__(self):

        x = ['%s Orbit State\n' % self.name]

        if self._r:
            x.append('r:     %0.4f km\n' % self._r)

        if self._r_rth is not None:
            x.append('r-rth:  %s km\n' % str(self._r_rth))

        if self._r_eph is not None:
            x.append('r-eph:  %s km\n' % str(self._r_eph))

        if self._r_xyz is not None:
            x.append('r-xyz:  %s km\n' % str(self._r_xyz))

        if self._v:
            x.append('v:     %0.4f km/s\n' % self._v)

        if self._v_rth is not None:
            x.append('v-rth:  %s km/s\n' % str(self._v_rth))

        if self._v_eph is not None:
            x.append('v-eph:  %s km/s\n' % str(self._v_eph))

        if self._v_xyz is not None:
            x.append('v-xyz:  %s km/s\n' % str(self._v_xyz))

        if self._ta is not None:
            x.append('ta:    %0.4f deg\n' % np.rad2deg(self._ta))
            x.append('ta:    %0.4f rad\n' % self._ta)

            if self._ta < 0:
                pos_ta = 2 * math.pi + self._ta
                x.append('ta:    %0.4f deg\n' % np.rad2deg(pos_ta))
                x.append('ta:    %0.4f rad\n' % pos_ta)

        if self._E is not None:
            x.append('E:     %0.4f deg\n' % np.rad2deg(self._E))
            x.append('E:     %0.4f rad\n' % self._E)

            if self._E < 0:
                pos_E = 2 * math.pi + self._E
                x.append('E:     %0.4f deg\n' % np.rad2deg(pos_E))
                x.append('E:     %0.4f rad\n' % pos_E)

        if self._M is not None:
            x.append('M:     %0.4f deg\n' % np.rad2deg(self._M))
            x.append('M:     %0.4f rad\n' % self._M)

        if self._fpa is not None:
            x.append('fpa:   %0.4f deg\n' % np.rad2deg(self._fpa))
            x.append('fpa:   %0.4f rad\n' % self._fpa)

        if self._theta is not None:
            x.append('theta: %0.4f deg\n' % np.rad2deg(self._theta))
            x.append('theta: %0.4f rad\n' % self._theta)

        if self._ttp is not None:
            x.append('ttp:   %0.4f s\n' % self._ttp)
            x.append('ttp:   %0.4f hr\n' % (self._ttp / 3600.0))

        if self._dcm_er is not None:
            x.append('C_er:  \n' + str(self._dcm_er) + '\n')

        if self._dcm_rv is not None:
            x.append('C_rv:  \n' + str(self._dcm_rv) + '\n')

        if self._dcm_ri is not None:
            x.append('C_ri:  \n' + str(self._dcm_ri) + '\n')

        return ''.join(x)


class ConicOrbit(object):

    def __init__(self, central_body, name=""):
        self.central_body = central_body
        self.name = name

        self._e = None
        self._a = None
        self._b = None
        self._rp = None
        self._ra = None
        self._p = None
        self._period = None
        self._se = None
        self._h = None
        self._h_xyz = None
        self._n = None
        self._asc_node = None
        self._arg_periapsis = None
        self._i = None
        self._vc = None

    def set_vars(self):
        self.e = None
        self.a = None
        self.b = None
        self.rp = None
        self.ra = None
        self.p = None
        self.period = None
        self.se = None
        self.h = None
        self.h_xyz = None
        self.n = None
        self.asc_node = None
        self.arg_periapsis = None
        self.i = None

    def set_circular(self):
        self.e = 0

        if self._a:
            self._p = self._a
            self._vc = np.sqrt(self.central_body.mu / self._a)

    @property
    def i(self):
        return self._i

    @i.setter
    def i(self, i=None):
        if self.i is not None:
            return

        if i is not None:
            self._i = i
        elif self.h_xyz is not None:
            self._i = np.arccos(self.h_xyz[2])
        else:
            return

        self.set_vars()

    @property
    def asc_node(self):
        return self._asc_node

    @asc_node.setter
    def asc_node(self, asc_node=None):
        if self._asc_node is not None:
            return

        if asc_node is not None:
            self._asc_node = asc_node
        elif self.h_xyz is not None and self.i is not None:
            sin_val = np.arcsin(self.h_xyz[0] / np.sin(self.i))
            cos_val = np.arccos(-self.h_xyz[1] / np.sin(self.i))
            self._asc_node = self._common_val(sin_val, cos_val)
        else:
            return

        self.set_vars()

    @property
    def arg_periapsis(self):
        return self._arg_periapsis

    @arg_periapsis.setter
    def arg_periapsis(self, arg_periapsis=None):
        if self._arg_periapsis is not None:
            return

        if arg_periapsis is not None:
            self._arg_periapsis = arg_periapsis
        else:
            return

        self.set_vars()

    @property
    def p(self):
        return self._p

    @p.setter
    def p(self, p=None):
        if self.p is not None:
            return

        if p:
            self._p = p

        elif self.a is not None and self.e is not None and self.is_ellipse():
            self._p = self.a * (1 - self.e ** 2)
        elif self.a is not None and self.e is not None and self.is_hyperbola():
            self._p = np.abs(self.a) * (self.e ** 2 - 1)
        elif self.h is not None:
            self._p = self.h ** 2 / self.central_body.mu
        else:
            return

        self.set_vars()

    @property
    def vc(self):
        return self._vc

    @vc.setter
    def vc(self, vc):
        self._vc = vc

    def get_rp(self):
        st = KeplarianState(self, "Periapsis")
        st.r = self.rp
        st.ta = 0.0
        st.set_ascending()
        return st

    def get_ra(self):
        st = KeplarianState(self, "Apoapsis")
        st.r = self.ra
        st.a = np.pi
        st.set_ascending()
        return st

    def is_circular(self):
        return self.e is not None and self._e == 0

    def is_ellipse(self):
        return (self.e is not None and self.e < 1.0) or (self._a is not None and self._a > 0)

    def is_hyperbola(self):
        return (self.e is not None and self.e > 1.0) or (self._a is not None and self._a < 0)

    def __repr__(self):
        x = ['%s Orbit Info\n' % self.name]

        if self._a is not None:
            x.append('a:      %0.4f km\n' % self._a)

        if self._b is not None:
            x.append('b:      %0.4f km\n' % self._b)

        if self.e is not None:
            x.append('e:      %0.4f\n' % self.e)

        if self._rp is not None:
            x.append('rp:     %0.4f km\n' % self._rp)

        if self._ra is not None:
            x.append('ra:     %0.4f km\n' % self._ra)

        if self._p is not None:
            x.append('p:      %0.4f km\n' % self._p)

        if self._period is not None:
            x.append('period: %0.4f s\n' \
                     'period: %0.4f hr\n' \
                     'period: %0.4f days\n' % (
                         self._period,
                         self._period / 3600,
                         self._period / 3600 / 24,
                     ))

        if self._h is not None:
            x.append('h:      %0.4f km^2/s\n' % self._h)

        if self._h_xyz is not None:
            x.append('h_xyz:  %s km^2/s\n' % str(self._h_xyz))

        if self._se is not None:
            x.append('se:     %0.4e km^2/s^2\n' % self._se)

        if self.vc is not None:
            x.append('vc:     %0.4f km/s\n' % self.vc)

        if self._i is not None:
            x.append('i:      %0.4f rad\n' % self._i)
            x.append('i:      %0.4f deg\n' % np.rad2deg(self._i))

        if self._asc_node is not None:
            x.append('Omega:  %0.4f rad\n' % self._asc_node)
            x.append('Omega:  %0.4f deg\n' % np.rad2deg(self._asc_node))

        if self._arg_periapsis is not None:
            x.append('omega:  %0.4f rad\n' % self._arg_periapsis)
            x.append('omega:  %0.4f deg\n' % np.rad2deg(self._arg_periapsis))

        return ''.join(x)


class LambertArc(object):
    def __init__(self, central_body, r1, r2, ta):
        self.r1 = r1
        self.r2 = r2

        self.ta = ta
        self.body = central_body

        self.c = conics_utils.loc(r1, r2, self.ta)
        self.s = 0.5 * (r1 + r2 + self.c)
        self.a_m = self.s / 2.
        self.alpha_m = 2 * np.arcsin(np.sqrt(self.s / (2 * self.a_m)))
        self.beta_m = 2 * np.arcsin(np.sqrt((self.s - self.c) / (2 * self.a_m)))
        self.TOF_m = np.sqrt(self.a_m ** 3 / self.body.mu) * \
                     (self.alpha_m - self.beta_m - np.sin(self.alpha_m) + np.sin(self.beta_m))
        self.p_m = (4. * self.a_m * (self.s - r1) * (self.s - r2)) / (self.c ** 2) * \
                   np.sin((self.alpha_m + self.beta_m) / 2.) ** 2

        self.transfer_orbit = ConicOrbit(self.body, 'Lambert Transfer')
        self.transfer_orbit.a = self.a_m
        self.transfer_orbit.p = self.p_m

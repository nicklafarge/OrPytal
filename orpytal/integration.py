import numpy as np
from scipy import integrate
from orpytal import OrbitType, units


def integrate_orbit(orbit, tol=1e-13, method="RK45"):
    if orbit.type() == OrbitType.Elliptic:
        start_state = orbit.get_state(ta=0)
        tend = orbit.period.m
    elif orbit.type() == OrbitType.Hyperbolic:
        start_state = orbit.get_state(ta=-orbit.ta_inf + 0.3)
        # start_state = orbit.get_state(ta=0)
        tend = -start_state.t_since_rp.m * 2
        # tend=((5*units.hr).to('s')).m
    else:
        raise ValueError("Only handles elliptic and hyperbolic")

    if orbit.angles_set():
        ic = np.concatenate((start_state.position.inertial(start_state).value.m,
                             start_state.velocity.inertial(start_state).value.m))
    else:
        ic = np.concatenate((start_state.position.perifocal(start_state).value.m,
                             start_state.velocity.perifocal(start_state).value.m))

    sol = integrate.solve_ivp(lambda t, y: two_body_eom(t, y, orbit.central_body),
                              (0, tend),
                              ic,
                              method=method,
                              dense_output=True,
                              # t_eval=np.linspace(0, orbit.period.m, 200),
                              rtol=tol,
                              atol=tol,
                              )

    return sol


def two_body_eom(t, y, body):
    r_vec = y[0:3]
    r = np.linalg.norm(r_vec)
    vel = y[3:6]
    accel = - (body.mu.m / r ** 3) * r_vec
    return np.concatenate((vel, accel))

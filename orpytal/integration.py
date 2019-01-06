import numpy as np
from scipy import integrate


def integrate_orbit(orbit, tol=1e-12, method="RK45"):
    periapsis_state = orbit.get_state(ta=0)

    if orbit.angles_set():
        ic = np.concatenate((periapsis_state.position.inertial().value.m, periapsis_state.velocity.inertial().value.m))
    else:
        ic = np.concatenate(
            (periapsis_state.position.orbit_fixed().value.m, periapsis_state.velocity.orbit_fixed().value.m))

    return integrate.solve_ivp(lambda t, y: two_body_eom(t, y, orbit.central_body),
                               (0, orbit.period.m),
                               ic,
                               method=method,
                               dense_output=True,
                               # t_eval=np.linspace(0, orbit.period.m, 200),
                               rtol=tol,
                               atol=tol,
                               )


def two_body_eom(t, y, body):
    r_vec = y[0:3]
    r = np.linalg.norm(r_vec)
    vel = y[3:6]
    accel = - (body.mu.m / r ** 3) * r_vec
    return np.concatenate((vel, accel))

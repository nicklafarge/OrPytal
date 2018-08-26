from orpytal.ksp import ksp_planet_constants
from orpytal import plotting
from matplotlib import pyplot as plt
from orpytal import matlab_plotting as mplt

kerbin = ksp_planet_constants.kerbin
mun = ksp_planet_constants.mun

plt.figure(1)
plotting.plot_orbit_inertial(ksp_planet_constants.moho.orbit)
plotting.plot_orbit_inertial(ksp_planet_constants.eve.orbit)
plotting.plot_orbit_inertial(ksp_planet_constants.kerbin.orbit)
plotting.plot_orbit_inertial(ksp_planet_constants.duna.orbit)
plotting.plot_orbit_inertial(ksp_planet_constants.dres.orbit)
plotting.plot_orbit_inertial(ksp_planet_constants.jool.orbit)
plotting.plot_orbit_inertial(ksp_planet_constants.eeloo.orbit)
plt.title('Kerbal Planets')
plt.show(block=False)

plotting.plot_orbit_inertial_3d(ksp_planet_constants.mun.orbit)
plotting.plot_orbit_inertial_3d(ksp_planet_constants.minimus.orbit)
mplt.title('Kerbin Moons')
mplt.plot_primary(ksp_planet_constants.kerbin)
mplt.legend()

# plt.figure(2)
# plotting.plot_orbit_inertial(ksp_planet_constants.laythe.orbit)
# plotting.plot_orbit_inertial(ksp_planet_constants.vall.orbit)
# plotting.plot_orbit_inertial(ksp_planet_constants.tylo.orbit)
# plotting.plot_orbit_inertial(ksp_planet_constants.bop.orbit)
# plotting.plot_orbit_inertial(ksp_planet_constants.pol.orbit)
# plt.title('Moons of Jool')
# plt.show(block=False)

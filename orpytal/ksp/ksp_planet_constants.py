##########################################################
# Reference: https://wiki.kerbalspaceprogram.com/wiki/Category:Celestials
##########################################################


########### Standard ###########

########### Local ###########
from common import units
from planet_constants import CentralBody

########### External ###########

##########################################################
# Values
##########################################################

G = 6.67408e-11 / (1000 ** 3) * units.km ** 3 / (units.kg * units.s ** 2)

##########################################################
# Planets
##########################################################

mu_units = (units.km ** 3 / units.s ** 2)

BODIES = {}

################## Sun ##################

BODIES['KERBOL'] = CentralBody(
    name='Kerbol',
    radius=26160000 * units.m,
    mu=1.1723328e18 * units('m^3/s^2')
)

################## Moho ##################

BODIES['MOHO'] = CentralBody(
    name='Moho',
    radius=250000 * units.m,
    mu=1.6860938e11 * units('m^3/s^2'),
    a=5263138304 * units.m,
    period=5657995 * units.s,
    e=0.2 * units.dimensionless,
    i=7 * units.deg,
    parent=BODIES['KERBOL'],
    arg_periapsis=15 * units.deg,
    ascending_node=70 * units.deg

)
################## Eve ##################

BODIES['EVE'] = CentralBody(
    name='Eve',
    radius=700000 * units.m,
    mu=8.1717302e12 * units('m^3/s^2'),
    a=9832684544 * units.m,
    period=5657995 * units.s,
    e=0.01 * units.dimensionless,
    i=2.1 * units.deg,
    parent=BODIES['KERBOL'],
    arg_periapsis=0 * units.deg,
    ascending_node=15 * units.deg
)

BODIES['GILLY'] = CentralBody(
    name='Gilly',
    radius=31500000 * units.m,
    mu=8289449.8 * units('m^3/s^2'),
    a=31500000 * units.m,
    period=388587 * units.s,
    e=0.55 * units.dimensionless,
    i=12 * units.degrees,
    parent=BODIES['EVE'],
    arg_periapsis=10 * units.deg,
    ascending_node=80 * units.deg
)

################## Kerbin ##################

BODIES['KERBIN'] = CentralBody(
    name='Kerbin',
    radius=600000 * units.m,
    mu=3.5316000e12 * units('m^3/s^2'),
    a=13599840256 * units.m,
    period=9203545 * units.s,
    e=0 * units.dimensionless,
    i=0 * units.deg,
    parent=BODIES['KERBOL'],
    arg_periapsis=0 * units.deg,
    ascending_node=0 * units.deg
)

BODIES['MUN'] = CentralBody(
    name='Mun',
    radius=200000 * units.m,
    mu=6.5138398e10 * units('m^3/s^2'),
    a=12000000 * units.m,
    period=138984 * units.s,
    e=0 * units.dimensionless,
    i=0 * units.degrees,
    parent=BODIES['KERBIN'],
    arg_periapsis=0 * units.deg,
    ascending_node=0 * units.deg
)

BODIES['MINIMUS'] = CentralBody(
    name='Minimus',
    radius=60000 * units.m,
    mu=1.7658000e9 * units('m^3/s^2'),
    a=47000000 * units.m,
    period=1077311 * units.s,
    e=0 * units.dimensionless,
    i=6 * units.degrees,
    parent=BODIES['KERBIN'],
    arg_periapsis=38 * units.deg,
    ascending_node=78 * units.deg
)

################## Duna ##################

BODIES['DUNA'] = CentralBody(
    name='Duna',
    radius=320000 * units.m,
    mu=3.0136321e11 * units('m^3/s^2'),
    a=20726155264 * units.m,
    period=17315400 * units.s,
    e=0.051 * units.dimensionless,
    i=0.06 * units.deg,
    parent=BODIES['KERBOL'],
    arg_periapsis=0 * units.deg,
    ascending_node=135.5 * units.deg
)

BODIES['IKE'] = CentralBody(
    name='Ike',
    radius=130000 * units.m,
    mu=1.8568369e10 * units('m^3/s^2'),
    a=3200000 * units.m,
    period=65518 * units.s,
    e=0.03 * units.dimensionless,
    i=0.2 * units.degrees,
    parent=BODIES['DUNA'],
    arg_periapsis=0 * units.deg,
    ascending_node=0 * units.deg
)

################## Dres ##################

BODIES['DRES'] = CentralBody(
    name='Dres',
    radius=138000 * units.m,
    mu=2.1484489e10 * units('m^3/s^2'),
    a=40839348203 * units.m,
    period=47893063 * units.s,
    e=0.145 * units.dimensionless,
    i=5 * units.deg,
    parent=BODIES['KERBOL'],
    arg_periapsis=90 * units.deg,
    ascending_node=280 * units.deg
)

################## Jool ##################

BODIES['JOOL'] = CentralBody(
    name='Jool',
    radius=6000000 * units.m,
    mu=2.8252800e14 * units('m^3/s^2'),
    a=68773560320 * units.m,
    period=104661432 * units.s,
    e=0.05 * units.dimensionless,
    i=1.304 * units.deg,
    parent=BODIES['KERBOL'],
    arg_periapsis=0 * units.deg,
    ascending_node=52 * units.deg
)

BODIES['LAYTHE'] = CentralBody(
    name='Laythe',
    radius=500000 * units.m,
    mu=1.9620000e12 * units('m^3/s^2'),
    a=27184000 * units.m,
    period=52981 * units.s,
    e=0 * units.dimensionless,
    i=0 * units.deg,
    parent=BODIES['JOOL'],
    arg_periapsis=0 * units.deg,
    ascending_node=0 * units.deg
)
BODIES['VALL'] = CentralBody(
    name='Vall',
    radius=300000 * units.m,
    mu=2.0748150e11 * units('m^3/s^2'),
    a=43152000 * units.m,
    period=105962 * units.s,
    e=0 * units.dimensionless,
    i=0 * units.deg,
    parent=BODIES['JOOL'],
    arg_periapsis=0 * units.deg,
    ascending_node=0 * units.deg
)
BODIES['TYLO'] = CentralBody(
    name='Jool',
    radius=600000 * units.m,
    mu=2.8252800e11 * units('m^3/s^2'),
    a=68500000 * units.m,
    period=211926 * units.s,
    e=0 * units.dimensionless,
    i=0.025 * units.deg,
    parent=BODIES['JOOL'],
    arg_periapsis=0 * units.deg,
    ascending_node=0 * units.deg
)
BODIES['BOP'] = CentralBody(
    name='Bop',
    radius=65000 * units.m,
    mu=2.4868349e9 * units('m^3/s^2'),
    a=128500000 * units.m,
    period=544507 * units.s,
    e=0.235 * units.dimensionless,
    i=15 * units.deg,
    parent=BODIES['JOOL'],
    arg_periapsis=25 * units.deg,
    ascending_node=10 * units.deg
)
BODIES['POL'] = CentralBody(
    name='Pol',
    radius=179890000 * units.m,
    mu=7.2170208e8 * units('m^3/s^2'),
    a=179890000 * units.m,
    period=901903 * units.s,
    e=0.171 * units.dimensionless,
    i=4.25 * units.deg,
    parent=BODIES['JOOL'],
    arg_periapsis=15 * units.deg,
    ascending_node=2 * units.deg
)

################## Eeloo ##################

BODIES['EELOO'] = CentralBody(
    name='Eeloo',
    radius=210000 * units.m,
    mu=7.4410815e10 * units('m^3/s^2'),
    a=90118820000 * units.m,
    period=156992048 * units.s,
    e=0.26 * units.dimensionless,
    i=6.15 * units.deg,
    parent=BODIES['KERBOL'],
    arg_periapsis=260 * units.deg,
    ascending_node=50 * units.deg
)

################## Easier Access ##################


kerbol = BODIES['KERBOL']

moho = BODIES['MOHO']

eve = BODIES['EVE']
gilly = BODIES['GILLY']

kerbin = BODIES['KERBIN']
mun = BODIES['MUN']
minimus = BODIES['MINIMUS']

duna = BODIES['DUNA']
ike = BODIES['IKE']

dres = BODIES['DRES']

jool = BODIES['JOOL']
laythe = BODIES['LAYTHE']
vall = BODIES['VALL']
tylo = BODIES['TYLO']
bop = BODIES['BOP']
pol = BODIES['POL']

eeloo = BODIES['EELOO']

########### Standard ###########
import os
import sys

########### Local ###########
from orpytal.common import units, utils
from orpytal import Orbit

########### External ###########
import untangle

##########################################################
# Values
##########################################################

AU = 149597870.7 * units.km
G = 6.67408e-11 / (1000 ** 3) * units.km ** 3 / (units.kg * units.s ** 2)

# km^3 / kg s^2


##########################################################
# Planets
##########################################################

mu_units = (units.km ** 3 / units.s ** 2)


class CentralBody(object):
    def __init__(self,
                 name,
                 radius,
                 mu,
                 a=None,
                 period=None,
                 e=None,
                 i=None,
                 parent=None,
                 id=None,
                 notes=None, **kwargs):
        self.name = name

        self.radius = radius * units.km if not isinstance(radius, units.Quantity) else radius.to(units.km)

        self.mu = mu * mu_units if not isinstance(mu, units.Quantity) else mu.to(mu_units)

        if a:
            self.a = a * units.km if not isinstance(a, units.Quantity) else a.to(units.km)

        if period:
            self.period = period * units.s if not isinstance(period, units.Quantity) else period.to(units.s)

        if e is not None:
            self.e = e * units.dimensionless if not isinstance(e, units.Quantity) else e.to(units.dimensionless)

        if i is not None:
            self.i = i * units.radians if not isinstance(i, units.Quantity) else i.to(units.radians)

        self.parent = parent

        if self.parent and isinstance(self.parent, CentralBody):
            self.orbit = Orbit(self.parent, '{}'.format(self.name),
                               a=self.a,
                               period=self.period,
                               e=self.e,
                               i=self.i,
                               **kwargs)
        if id:
            self.id = id

        if notes:
            self.notes = notes

    @classmethod
    def from_xml(cls, xml):
        id = int(xml.id.cdata)
        e = float(xml.ecc.cdata) * units.dimensionless
        mu = float(xml.gm.cdata) * (units.km ** 3 / units.s ** 2)
        a = float(xml.circ_r.cdata) * units.km
        i = float(xml.inc.cdata) * units.degree
        radius = float(xml.radius.cdata) * units.km
        name = xml.name.cdata
        notes = xml.notes.cdata
        parent = xml.parent.cdata

        return CentralBody(
            name,
            radius,
            mu,
            id=id,
            e=e,
            i=i,
            a=a,
            notes=notes,
            parent=None if parent == 'N/A' else parent
        )

    @classmethod
    def from_dict(cls, dict):
        return cls(name=dict['name'],
                   radius=dict['radius'],
                   mu=dict['mu'],
                   a=None if 'iad' not in dict else dict['a'],
                   period=None if 'period' not in dict else dict['period'],
                   e=None if 'e' not in dict else dict['e'],
                   i=None if 'i' not in dict else dict['i'],
                   parent=None if 'parent' not in dict else dict['parent'],
                   id=None if 'id' not in dict else dict['id'],
                   notes=None if 'notes' not in dict else dict['notes'],
                   )

    @classmethod
    def from_pickle(cls, filename):
        return cls.from_dict(utils.load_pickle(filename))

    def pickle(self, filename):
        utils.pickle_dict(filename, self.to_pickle_dict())

    def to_pickle_dict(self):
        dict = {
            'name': self.name,
            'radius': self.radius,
            'mu': self.mu
        }

        if self.a:
            dict['a'] = self.a
        if self.period:
            dict['period'] = self.period
        if self.e:
            dict['e'] = self.e
        if self.i:
            dict['i'] = self.i
        if self.parent:
            dict['parent'] = self.parent
        if self.id:
            dict['id'] = self.id
        if self.notes:
            dict['notes'] = self.notes

        return dict


class BodiesDict(object):
    def __init__(self, xml_file):

        self.__dict__ = {}

        all_bodies_xml = untangle.parse(xml_file)
        for body_xml in all_bodies_xml.body_data.body:
            body = CentralBody.from_xml(body_xml)
            self[body.name.lower()] = body

        for key, value in self.__dict__.items():
            if value.parent:
                value.parent = self[value.parent]

    def __getitem__(self, item):
        return self.__dict__[item.lower()]

    def __setitem__(self, key, value):
        self.__dict__[key.lower()] = value

    def __repr__(self):
        return repr(self.__dict__)

    def __len__(self):
        return len(self.__dict__)

    def __delitem__(self, key):
        del self.__dict__[key]

    def clear(self):
        return self.__dict__.clear()

    def copy(self):
        return self.__dict__.copy()

    def has_key(self, k):
        return k in self.__dict__

    def update(self, *args, **kwargs):
        return self.__dict__.update(*args, **kwargs)

    def keys(self):
        return self.__dict__.keys()

    def values(self):
        return self.__dict__.values()

    def items(self):
        return self.__dict__.items()

    def pop(self, *args):
        return self.__dict__.pop(*args)

    def __cmp__(self, dict_):
        return self.__cmp__(self.__dict__, dict_)

    def __contains__(self, item):
        return item in self.__dict__

    def __iter__(self):
        return iter(self.__dict__)

BODIES = BodiesDict(os.path.join(os.path.dirname(__file__), 'data/body_data.xml'))

# Add periods
BODIES['MERCURY'].period = 7600531.965 * units.seconds
BODIES['VENUS'].period = 19414144.04 * units.seconds
BODIES['MOON'].period = 236059 * units.seconds
BODIES['EARTH'].period = 31574658.69 * units.seconds
BODIES['PHOBOS'].period = 27552.96 * units.seconds
BODIES['DEIMOS'].period = 109071.36 * units.seconds
BODIES['MARS'].period = 59359794.38 * units.seconds
BODIES['IO'].period = 152853.5047 * units.seconds
BODIES['EUROPA'].period = 306822.0384 * units.seconds
BODIES['GANYMEDE'].period = 618153.3757 * units.seconds
BODIES['JUPITER'].period = 374869015.7 * units.seconds
BODIES['TITAN'].period = 1377648 * units.seconds
BODIES['SATURN'].period = 929260480.9 * units.seconds
BODIES['TITANIA'].period = 752218.6176 * units.seconds
BODIES['URANUS'].period = 2656887761 * units.seconds
BODIES['TRITON'].period = .5877 * 3600 * units.seconds
BODIES['NEPTUNE'].period = 5193323332 * units.seconds
BODIES['CHARON'].period = 551856.7066 * units.seconds
BODIES['PLUTO'].period = 7753066072 * units.seconds

BODIES['CERES'] = CentralBody(
    name='Ceres',
    radius=476.20,
    mu=63.2,
    a=413690604,
    period=145238090.6,
    e=0.0757972598,
    i=10.593981 * units.degrees,
    parent=BODIES['SUN']
)

BODIES['CALLISTO'] = CentralBody(
    name='Callisto',
    radius=2410.30,
    mu=7179.289,
    a=1882700,
    period=1441931.19,
    e=0.0074,
    i=0.192 * units.degrees,
    parent=BODIES['JUPITER']
)

BODIES_532 = {}

################## Sun ##################

BODIES_532['SUN'] = CentralBody(
    name='Sun',
    radius=695990.00,
    mu=132712200000.00
)

################## Mercury ##################

BODIES_532['MERCURY'] = CentralBody(
    name='Mercury',
    radius=2439.70,
    mu=22032.09,
    a=57909043,
    period=7600531.965,
    e=0.205645,
    i=7.004161 * units.degrees

)

################## Venus ##################

BODIES_532['VENUS'] = CentralBody(
    name='Venus',
    radius=6051.80,
    mu=324858.592079,
    a=108208877,
    period=19414144.04,
    e=0.006766272,
    i=3.394656 * units.degrees
)
################## Earth ##################


BODIES_532['EARTH'] = CentralBody(
    name='Earth',
    radius=6378.137,
    mu=398600.4418,
    a=149649952,
    period=31574658.69,
    e=0.01621148,
    i=0.0008611136 * units.degrees
)

BODIES_532['MOON'] = CentralBody(
    name='Moon',
    radius=1737.40,
    mu=4902.801076,
    a=384400,
    period=2360592,
    e=0.0543866,
    i=5.16000 * units.degrees,
    parent=BODIES_532['EARTH']
)

################## Mars ##################

BODIES_532['MARS'] = CentralBody(
    name='Mars',
    radius=3396.19,
    mu=42828.3719012840,
    a=227953016,
    period=59359794.38,
    e=0.0933555537,
    i=1.848489
)

BODIES_532['PHOBOS'] = CentralBody(
    name='Phobos',
    radius=11.10,
    mu=0.0007112,
    a=9376,
    period=27552.96,
    e=0.01510000,
    i=1.075000,
    parent=BODIES_532['MARS']
)

BODIES_532['DEIMOS'] = CentralBody(
    name='Deimos',
    radius=6.20,
    mu=0.0000985,
    a=23458,
    period=109071.36,
    e=0.0002,
    i=1.78800,
    parent=BODIES_532['MARS']
)

################## Jupiter ##################

BODIES_532['JUPITER'] = CentralBody(
    name='Jupiter',
    radius=71492.00,
    mu=126686535,
    a=779067093,
    period=374869015.7,
    e=0.0490681913,
    i=1.303603 * units.degrees
)

BODIES_532['IO'] = CentralBody(
    name='Io',
    radius=1821.60,
    mu=5959.916,
    a=421800,
    period=152853.5047,
    e=0.0041,
    i=0.036 * units.degrees,
    parent=BODIES_532['JUPITER']
)

BODIES_532['EUROPA'] = CentralBody(
    name='Europa',
    radius=1560.80,
    mu=3202.739,
    a=671100,
    period=306822.0384,
    e=0.0094,
    i=0.466 * units.degrees,
    parent=BODIES_532['JUPITER']
)

BODIES_532['GANYMEDE'] = CentralBody(
    name='Ganymede',
    radius=2631.20,
    mu=9887.834,
    a=1070400,
    period=618153.3757,
    e=0.0013,
    i=0.177 * units.degrees,
    parent=BODIES_532['JUPITER']
)

BODIES_532['CALLISTO'] = CentralBody(
    name='Callisto',
    radius=2410.30,
    mu=7179.289,
    a=1882700,
    period=1441931.19,
    e=0.0074,
    i=0.192 * units.degrees,
    parent=BODIES_532['JUPITER']
)

################## Saturn ##################

BODIES_532['SATURN'] = CentralBody(
    name='Saturn',
    radius=60268.00,
    mu=37931284,
    a=1426647631,
    period=929260480.9,
    e=0.05483734,
    i=2.488476 * units.degrees
)

BODIES_532['TITAN'] = CentralBody(
    name='Titan',
    radius=2574.730,
    mu=8978.1382,
    a=1221865,
    period=1377648,
    e=0.0288,
    i=0.312 * units.degrees,
    parent=BODIES_532['SATURN']
)

################## Uranus ##################

BODIES_532['URANUS'] = CentralBody(
    name='Uranus',
    radius=25559.00,
    mu=5793965.66393928,
    a=2873682806,
    period=2656887761,
    e=0.046796657,
    i=0.77242 * units.degrees
)

BODIES_532['TITANIA'] = CentralBody(
    name='Titania',
    radius=788.90,
    mu=228.2,
    a=436300,
    period=752218.6176,
    e=0.0011,
    i=0.079 * units.degrees,
    parent=BODIES_532['URANUS']
)

################## Neptune ##################


BODIES_532['NEPTUNE'] = CentralBody(
    name='Neptune',
    radius=24764,
    mu=6835110.00,
    a=4492499814,
    period=5193323332,
    e=0.010016637,
    i=1.767227 * units.degrees
)

BODIES_532['TRITON'] = CentralBody(
    name='Triton',
    radius=1350.00,
    mu=1427.598,
    a=354759,
    period=.5877 * 3600,
    e=0.000016,
    i=156.86500 * units.degrees,
    parent=BODIES_532['NEPTUNE']
)

################## Pluto ##################


BODIES_532['PLUTO'] = CentralBody(
    name='Pluto',
    radius=1195.00,
    mu=873.7674473980320,
    a=5868125070,
    period=7753066072,
    e=0.247576314,
    i=17.325236 * units.degrees
)

BODIES_532['CHARON'] = CentralBody(
    name='Charon',
    radius=603.60,
    mu=102.3,
    a=17536,
    period=551856.7066,
    e=0.0022,
    i=0.001 * units.degrees,
    parent=BODIES_532['PLUTO']
)

################## Ceres ##################

BODIES_532['CERES'] = CentralBody(
    name='Ceres',
    radius=476.20,
    mu=63.2,
    a=413690604,
    period=145238090.6,
    e=0.0757972598,
    i=10.593981 * units.degrees

)

for k, b in BODIES.items():
    setattr(sys.modules[__name__], b.name.lower(), b)

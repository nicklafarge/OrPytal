import orpytal
from orpytal import Orbit, KeplarianState, units

DEFAULT_LENGTH = 17
DEFAULT_PRECISION = 10
DEFAULT_PRECISION_TOTAL_LENGTH = DEFAULT_PRECISION + 7
DEFAULT_UNIT_LENGTH = 10

DEFAULT_VALUE = '.'
DEFAULT_FILLER = DEFAULT_VALUE * DEFAULT_LENGTH
DEFAULT_UNIT_FILLER = DEFAULT_VALUE * DEFAULT_UNIT_LENGTH

LABEL_FMT = '{{0: >{}}}'.format(DEFAULT_LENGTH)
VALUE_FMT = '{{:{}.{}e}}'.format(DEFAULT_PRECISION_TOTAL_LENGTH, DEFAULT_PRECISION)
UNIT_FMT = '{{0: >{}}}'.format(DEFAULT_UNIT_LENGTH)

VALUE_SEPARATOR_LENGTH = 5
VALUE_SEPARATOR = " " * VALUE_SEPARATOR_LENGTH

ELEMENT_LENGTH = DEFAULT_LENGTH + DEFAULT_PRECISION_TOTAL_LENGTH + DEFAULT_UNIT_LENGTH + 5

LINE_WIDTH = ELEMENT_LENGTH * 2 + 1 * VALUE_SEPARATOR_LENGTH
H_DIVIDER_CHARACTER = "-"
SECTION_DIVIDER = H_DIVIDER_CHARACTER * LINE_WIDTH


def convert_unit_str(unit_str):
    # General
    unit_str = unit_str.replace(" ", "")
    unit_str = unit_str.replace("**", "^")

    # Specific
    unit_str = unit_str.replace("kilometer", "km")
    unit_str = unit_str.replace("second", "s")
    unit_str = unit_str.replace("dimensionless", "nd")
    unit_str = unit_str.replace("radian", "rad")
    unit_str = unit_str.replace("degree", "deg")

    return unit_str


def create_key_value(label, value):
    val = DEFAULT_FILLER if not value else VALUE_FMT.format(value.m)
    u = DEFAULT_UNIT_FILLER if not value else UNIT_FMT.format(convert_unit_str(str(value.u)))
    return "{}: {} {} |".format(LABEL_FMT.format(label), val, u)


def create_key_value_line(labels, values):
    output = []
    for l, v in zip(labels, values):
        output.append(create_key_value(l, v))
        output.append(VALUE_SEPARATOR)
    return "".join(output)


def create_key_value_line_params(*orbit_params):
    labels = []
    values = []
    for p in orbit_params:
        labels.append(p.name)
        values.append(p.value)

    return create_key_value_line(labels, values)


def output_orbit(orbit):
    output = []

    # Construct Header
    header = []
    header.append(SECTION_DIVIDER)
    header.append("Orbit Data".center(LINE_WIDTH, " "))
    header.append(SECTION_DIVIDER)
    header.append("Meta Information".center(LINE_WIDTH, " "))
    header.append("{:>12}: {}".format("Central Body", orbit.central_body.name))
    if orbit.name:
        header.append("{:>12}: {}".format("Orbit Name", orbit.name))

    header.append(SECTION_DIVIDER)

    output.append("\n".join(header))

    # Keplerian Elements
    output.append("Keplerian Elements".center(LINE_WIDTH, " ") + "\n")
    output.append(create_key_value_line_params(orbit._e, orbit._i))
    output.append(create_key_value_line_params(orbit._a, orbit._ascending_node))
    output.append(create_key_value_line_params(orbit._arg_periapsis))

    # Other Orbit Elements
    output.append("\n" + "Orbital Distances".center(LINE_WIDTH, " ") + "\n")
    output.append(create_key_value_line_params(orbit._rp, orbit._ra))
    output.append(create_key_value_line_params(orbit._p, orbit._b))

    output.append("\n" + "Orbit Rate and Period".center(LINE_WIDTH, " ") + "\n")
    output.append(create_key_value_line_params(orbit._n, orbit._period))

    output.append("\n" + "Other Orbital Elements".center(LINE_WIDTH, " ") + "\n")
    output.append(create_key_value_line_params(orbit._se, orbit._h))
    output.append(SECTION_DIVIDER)

    output_str = "\n".join(output)
    print(output_str)
    pass


def output_state(state):
    pass


if __name__ == "__main__":
    orbit = Orbit(orpytal.bodies.earth, name="My Dope Orbit", a=25000 * units.km, e=0.4, ascending_node=45 * units.deg,
                  inclination=16 * units.deg)
    output_orbit(orbit)
    pass

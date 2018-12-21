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
    val = DEFAULT_FILLER if value is None else VALUE_FMT.format(value.m)
    u = DEFAULT_UNIT_FILLER if value is None else UNIT_FMT.format(convert_unit_str(str(value.u)))
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


def orbit_parameters_output(orbit):
    output = []
    output.append("Orbit Parameters".center(LINE_WIDTH, " ") + "\n")
    output.append(create_key_value_line_params(orbit._e, orbit._i))
    output.append(create_key_value_line_params(orbit._a, orbit._raan))
    output.append(create_key_value_line_params(orbit._b, orbit._arg_periapsis))
    output.append(create_key_value_line_params(orbit._rp, orbit._ra))
    output.append(create_key_value_line_params(orbit._p, orbit._h))
    output.append(create_key_value_line_params(orbit._n, orbit._period))
    output.append(create_key_value_line_params(orbit._se))
    return output


def output_orbit(orbit):
    output = []
    output.append(SECTION_DIVIDER)

    # Construct Header
    header = []
    header.append(SECTION_DIVIDER)
    header.append("Orbit Data".center(LINE_WIDTH, " "))
    header.append(SECTION_DIVIDER)
    header.append("Meta Information".center(LINE_WIDTH, " "))
    header.append("{:>12}: {}".format("Orbit Name", DEFAULT_FILLER if not orbit.name else orbit.name))
    header.append("{:>12}: {}".format("Central Body", orbit.central_body.name))
    header.append("{:>12}: {}".format("Circular", DEFAULT_FILLER if orbit.e is None else orbit.circular()))
    header.append("{:>12}: {}".format("Equitorial", DEFAULT_FILLER if orbit.i is None else orbit.equitorial()))

    header.append(SECTION_DIVIDER)

    output.append("\n".join(header))
    output.append("\n".join(orbit_parameters_output(orbit)))
    output.append(SECTION_DIVIDER)
    output.append(SECTION_DIVIDER)

    output_str = "\n".join(output)

    return output_str


def output_state(state):
    orbit = state.orbit

    output = []
    output.append(SECTION_DIVIDER)

    # Construct Header
    header = []
    header.append(SECTION_DIVIDER)
    header.append("Keplarian State Data".center(LINE_WIDTH, " "))
    header.append(SECTION_DIVIDER)
    header.append("Meta Information".center(LINE_WIDTH, " "))
    header.append("{:>12}: {}".format("Orbit Name", DEFAULT_FILLER if not orbit.name else orbit.name))
    header.append("{:>12}: {}".format("State Name", DEFAULT_FILLER if not state.name else state.name))
    header.append("{:>12}: {}".format("Central Body", orbit.central_body.name))
    header.append("{:>12}: {}".format("Circular", DEFAULT_FILLER if orbit.e is None else orbit.circular()))
    header.append("{:>12}: {}".format("Equitorial", DEFAULT_FILLER if orbit.i is None else orbit.equitorial()))
    header.append("{:>12}: {}".format("Ascending", DEFAULT_FILLER if state._is_ascending is None else state._is_ascending))

    header.append(SECTION_DIVIDER)

    output.append("\n".join(header))

    # Keplarian State Data

    output.append("State Parameters".center(LINE_WIDTH, " ") + "\n")
    output.append(create_key_value_line_params(state._r, state._ta))
    output.append(create_key_value_line_params(state._v, state._E))
    output.append(create_key_value_line_params(state._fpa, state._M))
    output.append(create_key_value_line_params(state._t_since_rp, state._arg_latitude))
    output.append(SECTION_DIVIDER)
    output.append(SECTION_DIVIDER)

    # Orbit

    output.append("\n".join(orbit_parameters_output(orbit)))
    output.append(SECTION_DIVIDER)


    output_str = "\n".join(output)

    return output_str

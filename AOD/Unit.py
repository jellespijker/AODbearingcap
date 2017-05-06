from pint import UnitRegistry, set_application_registry

ureg = UnitRegistry(autoconvert_offset_to_baseunit=True)  # Allows for unit save calculations
Q_ = ureg.Quantity  # Allows for custom quantities to be registered
g = 9.80665 * ureg['m/s**2']  # Gravitational constant
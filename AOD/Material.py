from math import sqrt, sin, exp, tan, pi, cos

from AOD.Unit import *
import sys


class Material(object):
    """Material class, the parent of all used materials"""
    _rho = 1000. * ureg['kg/m**3']  # private density
    _T = 15 * ureg['degC']  # private temperature
    _rho_func = None  # private function for materials where density is a function of temperature
    _depth = 1. * ureg['m']

    def __str__(self):
        return r'Properties at ' + str(self.T.to('degC')) + '\n' + 'density: ' + str(self.rho)

    def __init__(self, rho=None, T=None, rho_func=None):
        """Constructor for Materials"""
        if rho is not None:
            self._rho = rho
        if T is not None:
            self._T = T
        else:
            self.T = 15
        if rho_func is not None:
            self._rho_func = rho_func

    @property
    def depth(self):
        return self._depth

    @depth.setter
    def depth(self, value):
        self._depth = value

    @property
    def T(self):
        """ Returns the current temperature of model"""
        return self._T

    @T.setter
    def T(self, value):
        """ Sets the new temperature, in degrees Celsius """
        if type(value) is type(ureg['degC']):
            self._T = value
        else:
            self._T = value * ureg['degC']

    @property
    def rho(self):
        """ Gets the density of a material"""
        return self._rho

    @rho.setter
    def rho(self, value):
        """ Sets the density of a material"""
        self._rho = value


class Fluid(Material):
    _P = 0. * ureg['Pa']  # private pressure of the fluid
    _mu = 0. * ureg['Pa*s']  # private dynamic viscosity
    _mu_func = None  # private dynamic viscosity

    def __init__(self, rho=None, rho_func=None, T=None, P=None, mu=None, mu_func=None):
        Material.__init__(self, T=T, rho=rho, rho_func=rho_func)
        if P is not None:
            self._P = P
        if mu is not None:
            self._mu = mu
        if mu_func is not None:
            self._mu_func = mu_func

    @property
    def P(self):
        return self._P

    @property
    def mu(self):
        """Returns dynamic viscosity, when mu function is specified the function 
        is applied, otherwise the private variable is returned"""
        if self._mu_func is not None:
            return self._mu_func(self.T.to('degC').magnitude)
        else:
            return self._mu

    @mu.setter
    def mu(self, value):
        """ Sets the dynamic viscosity"""
        self._mu = value

    @property
    def nu(self):
        """ returns the kinematic viscosity"""
        return (self.mu / self.rho).to('m**2/s')


class Air(Fluid):
    _R = 0. * ureg['J/(kg*K)']  # private specific gas constant for dry air

    def __init__(self, T=None, P=None, R=None):
        if T is None:
            T = 15. * ureg['degC']
        if P is None:
            P = 101.325e3 * ureg['Pa']
        if R is None:
            self._R = 287.05 * ureg['J/(kg*K)']
        Fluid.__init__(self, rho=1.225 * ureg['kg/m**3'],
                       rho_func=lambda p, r, t: p.to_base_units() / (r.to_base_units() * t.to('K')),
                       T=T, P=P, mu=1.789e-5 * ureg['Pa*s'])
        self.depth = float('inf') * ureg['m']

    @property
    def R(self):
        return self._R

    @R.setter
    def R(self, value):
        self._R = value

    @property
    def rho(self):
        """ Returns the density of a air, when the material has a
        rho_function specified, the density is calculated using the ideal gas law, otherwise it is
        taken from the private variable """
        if self._rho_func is not None:
            return self._rho_func(self.P, self.R, self.T)
        else:
            return self._rho


class Water(Fluid):
    """ Water material """

    def __init__(self, T=None, P=None):
        """ Constructor of the water class"""
        if T is not None:
            self._T = T
        self._rho = 999.7 * ureg['kg/m**3']
        # Water density dependency on temperature (validity 5<T_f<100C) (Matousek, 2004)
        self._rho_func = \
            lambda T: (999.7 - 0.10512 * (T - 10) - 0.005121 * (T - 10) ** 2 + 0.00001329 * (T - 10) ** 3) * ureg[
                'kg/m**3']
        # Water dynamic & kinematic viscosity as function of temperature (Matousek, 2004)
        self._mu_func = \
            lambda T: (0.10 / (2.1482 * ((T - 8.435) + sqrt(8078.4 + (T - 8.435) ** 2)) - 120)) * ureg['Pa*s']
        self.depth = 30. * ureg['m']

    @property
    def rho(self):
        """ Returns the density of water, when the material has a
        rho_function specified, the density is calculated, otherwise it is
        taken from the private variable """
        if self._rho_func is not None:
            return self._rho_func(self._T.to('degC').magnitude)
        else:
            return self._rho


class Soil(Material):
    """ Definitions for the different types of soils """
    _c = 3.e3 * ureg['Pa']  # private cohesion
    _phi = 0. * ureg['degree']  # private internal friction angle
    _k0 = 0.54 * ureg['dimensionless']  # private coeff of lateral earth press
    _delta = 0. * ureg['rad']  # private external friction angle
    _rho_ins = 1300. * ureg['kg/m**3']  # private In-situ bottom density
    _alpha = {0.2: 0.95,
              0.3: 0.77,
              0.4: 0.7,
              0.5: 0.65,
              0.6: 0.62,
              0.7: 0.6,
              0.8: 0.56,
              0.9: 0.55,
              1.: 0.53,
              1.1: 0.52,
              1.2: 0.5,
              1.3: 0.49,
              1.4: 0.48,
              1.5: 0.47,
              1.6: 0.42,
              1.7: 0.41,
              1.8: 0.41,
              1.9: 0.42,
              2.0: 0.41,
              2.1: 0.41,
              2.2: 0.4,
              2.3: 0.4,
              2.4: 0.4,
              2.5: 0.4,
              3.0: 0.39,
              4.0: 0.39}

    def __str__(self):
        return Material.__str__(self) + ' \nCohesion: ' + str(self.c) + ' \nint. friction angle (phi): ' + str(
            self.phi) + '\ncoeff of lateral earth pressure: ' + str(self.k0) + '\nexternal friction angle: ' + str(
            self.delta) + '\nIn-situ density: ' + str(self.rho_ins)

    def __init__(self, rho=None, T=None, rho_ins=None, c=None, phi=None, k0=None, delta=None):
        """ Constructor for soil materials"""
        Material.__init__(self, rho, T)
        if rho_ins is not None:
            self.rho_ins = rho_ins
        if c is not None:
            self.c = c
        if phi is not None:
            self.phi = phi
        if k0 is not None:
            self.k0 = k0
        if delta is not None:
            self.delta = delta
        self.depth = -float('inf') * ureg['m']

    def alpha(self, weight, area):
        ratio = (weight / area).item(0)
        prev_key = sys.float_info.min
        if ratio in self._alpha.keys():
            return self._alpha[ratio]
        else:
            if len(self._alpha) == 1 or list(self._alpha.keys())[0] > ratio:
                return list(self._alpha.values())[0]
            for key in self._alpha.items():
                if key[0] > ratio:
                    dR = ratio - prev_key
                    dAlpha = self._alpha[key[0]] - self._alpha[prev_key]
                    return self._alpha[prev_key] + dR * dAlpha / (key[0] - prev_key)
                else:
                    prev_key = key[0]
        return self._alpha[prev_key]

    def alpha_bulldozer(self, sigma_sb, sigma_cb):
        ratio = sigma_sb / sigma_cb
        alpha = (-0.0262 * ratio ** 3 + 0.223 * ratio ** 2 - 0.6143 * ratio + 0.9509) + 0.1
        return alpha

    def gamma(self, layers):
        """ Submerged soil weight """
        return g * (layers['Soil'].rho_ins - layers['Fluid'].rho)

    def sigma_h(self, load):
        """ Max Horizontal total soil stress assuming soil is at rest """
        return load * self.k0

    @property
    def c(self):
        """ Gets the soil cohesion """
        return self._c

    @c.setter
    def c(self, value):
        """ Sets the soil cohesion """
        self._c = value

    @property
    def phi(self):
        """ Gets the soil internal friction angle """
        return self._phi

    @phi.setter
    def phi(self, value):
        """ Sets the soil internal friction angle """
        self._phi = value.to('rad')

    @property
    def k0(self):
        """ Gets the soil coeff. of lateral earth pressure """
        return self._k0

    @k0.setter
    def k0(self, value):
        """ Sets the soil coeff. of lateral earth pressure"""
        self._k0 = value

    @property
    def delta(self):
        """ Gets the soil external friction angle """
        return self._delta

    @delta.setter
    def delta(self, value):
        """ Sets the soil external friction angle """
        self._delta = value.to('rad')

    @property
    def rho_ins(self):
        """ Gets the soil In-situ density """
        return self._rho_ins

    @rho_ins.setter
    def rho_ins(self, value):
        """ Sets teh soil In-situ density """
        self._rho_ins = value.to('kg/m**3')

    @property
    def N_q(self):
        """ Dimensionless constants in the Brinch-Hansen model (verruijt, 2009) """
        return (1 + sin(self.phi)) / (1 - sin(self.phi)) * exp(pi * tan(self.phi))

    @property
    def N_gamma(self):
        return 2 * (self.N_q - 1) * tan(self.phi)

    @property
    def N_c(self):
        if self.phi == 0.:
            return 2 * pi
        else:
            return (self.N_q - 1) * (1 / tan(self.phi))  # TODO check if cot is 1/tan(alpha)

    def S_c(self, B, L):
        return 1. + 0.2 * (B / L)

    def S_q(self, B, L):
        return 1. + (B / L) * sin(self.phi)

    def S_gamma(self, B, L):
        return 1. - 0.3 * (B / L)

    def i_c(self, p=None, t=None):
        """ Inclination factor currently not used """
        if p is not None and t is not None:
            ic = 1 - (t / (self.c * p * tan(self.phi)))
        return 1.

    def i_q(self, p=None, t=None):
        """ Inclination factor currently not used """
        iq = self.i_c(p, t) ** 2
        return 1.

    def i_gamma(self, p=None, t=None):
        """ Inclination factor currently not used """
        igamma = self.i_c(p, t) ** 3
        return 1.

    def p_allow(self, q, layers, B, L):
        """ Allowed load according to Brinch Hansen  """
        return self.i_c() * self.S_c(B, L) * self.c * self.N_c \
               + self.i_q() * self.S_q(B, L) * q * self.N_q \
               + self.i_gamma() * self.S_gamma(B, L) * 0.5 * self.gamma(layers) * B * self.N_gamma


class Silt(Soil):
    """ Predefined type of Soil, namely silt"""

    def __init__(self):
        """ Constructor """
        Soil.__init__(self, rho=2650. * ureg['kg/m**3'], rho_ins=1300. * ureg['kg/m**3'], c=3.e3 * ureg['Pa'],
                      phi=0. * ureg['degree'],
                      k0=0.54 * ureg['dimensionless'], delta=0. * ureg['degree'])


class Loose_clay(Soil):
    """ Predefined type of Soil, namely Loose clay"""

    def __init__(self):
        """ Constructor """
        Soil.__init__(self, rho=2650. * ureg['kg/m**3'], rho_ins=1400. * ureg['kg/m**3'], c=5.e3 * ureg['Pa'],
                      phi=0. * ureg['degree'],
                      k0=0. * ureg['dimensionless'], delta=0. * ureg['degree'])


class Packed_clay(Soil):
    """ Predefined type of Soil, namely Packed clay"""

    def __init__(self):
        """ Constructor """
        Soil.__init__(self, rho=2650. * ureg['kg/m**3'], rho_ins=1800. * ureg['kg/m**3'], c=10.e3 * ureg['Pa'],
                      phi=0. * ureg['degree'],
                      k0=1. * ureg['dimensionless'], delta=0. * ureg['degree'])


class River_clay(Soil):
    """ Predefined type of Soil, namely River clay used during the test with project 64120-R02"""

    def __init__(self):
        """ Constructor """
        Soil.__init__(self, rho=2650. * ureg['kg/m**3'], rho_ins=(1.86 * ureg['kg/l'] + 1.88 * ureg['kg/l']) / 2,
                      c=(4.2 * ureg['kPa'] + 3.1 * ureg['kPa']) / 2,
                      phi=0. * ureg['degree'],
                      k0=1. * ureg['dimensionless'], delta=0. * ureg['degree'])


class Sand(Soil):
    """ Predefined type of Soil, namely Packed clay"""

    def __init__(self):
        """ Constructor """
        Soil.__init__(self, rho=2650. * ureg['kg/m**3'], rho_ins=2000. * ureg['kg/m**3'], c=0. * ureg['Pa'],
                      phi=35. * ureg['degree'],
                      k0=0. * ureg['dimensionless'], delta=35. / 3. * ureg['degree'])

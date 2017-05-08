import numpy as np
from AOD.Material import *
from AOD.Unit import *

from AOD.Bot import *


class World(object):
    """ The representation of fluid and soil layer(s), and their interaction,
     currently only one soil layer is supported """
    Layers = {}

    _T = 15 * ureg['degC']  # private temperature of both fluid and soil

    def __init__(self, layers=None, T=15 * ureg['degC'], depths=None):
        """Constructor specifying the materials for layers, the model temperature and 
        the depth for each layer, as a dict. The soilbed, demarcation between soil and fluid
        is 0. [m], where downwards is specified as negative and upwards as positive"""
        if layers is None:
            self.Layers['Air'] = Air()
            self.Layers['Fluid'] = Water()
            self.Layers['Soil'] = Silt()
        self.T = T
        if depths is not None:
            self.layerdepths = depths

    @property
    def layerdepths(self):
        """ Getter for the depth of each layer, as a dict. The soilbed, demarcation between soil and fluid
        is 0. [m], where downwards is specified as negative and upwards as positive"""
        d = {}
        for key, l in self.Layers.items():
            d[key] = l.depth
        return d

    @layerdepths.setter
    def layerdepths(self, value):
        """Setter for the depth for each layer, as a dict. The soilbed, demarcation between soil and fluid
        is 0. [m], where downwards is specified as negative and upwards as positive"""
        for key, l in self.Layers.items():
            if key in value:
                l.depth = value[key]
            else:
                print(key + ' depth not passed to layers, assuming: ' + str(l.depth))

    @property
    def n(self):
        """ Gets the porosity of the soilbed """
        return (self.Layers['Soil'].rho - self.Layers['Soil'].rho_ins) / (
            self.Layers['Soil'].rho - self.Layers['Fluid'].rho)

    @property
    def T(self):
        """ Gets the temperature of the model """
        return self._T

    @T.setter
    def T(self, value):
        """ Sets the temperature of the model, and subsequently, that of the fluid
        and soil"""
        for key, l in self.Layers.items():
            l.T = value
        self._T = value

    @property
    def S(self):
        """ Gets the specific gravity / rel. density of the soilbed """
        return self.Layers['Soil'].rho / self.Layers['Fluid'].rho

    @property
    def gamma(self):
        """ Gets the submerged soil weight of the soilbed """
        return (g * (self.Layers['Soil'].rho_ins - self.Layers['Fluid'].rho)).to('N/m**3')


class Model(object):
    """ The complete model, with setup and solver"""
    world = World()
    bot = Bot()

    def __init__(self, world=None, bot=None):
        if world is not None:
            self.world = world
        if bot is not None:
            self.bot = bot

    def solve_sinkdepth(self, depth=None, resolution=None):
        max_sink_depth = 10. * ureg['m']
        if resolution is None:
            resolution = 1.e-3 * ureg['m']
        if depth is None:
            depth = np.arange(start=0., stop=max_sink_depth.magnitude, step=resolution.magnitude) * ureg['m']
        self.bot.depth = depth
        [B_acc, d_acc] = self.bot.Screw.cylinder.B_acc(depth=depth)
        force = self.bot.F(forcetype=ForceType.Wet, layers=self.world.Layers)
        force /= self.bot.no_screw
        p_load = force / (B_acc * self.bot.Screw.cylinder.l)
        gamma = self.world.Layers['Soil'].gamma(self.world.Layers)
        q = gamma * depth
        p_allow = self.world.Layers['Soil'].p_allow(q, self.world.Layers, B_acc, self.bot.Screw.cylinder.l)

        p_eps = np.sign(p_allow - p_load)
        sink_depth = -1. * ureg['m']
        load = -1. * ureg['Pa']
        for i in range(2, len(depth)):
            if p_eps[i] * p_eps[i - 1] == -1:
                sink_depth = round(depth[i], 3)
                load = p_load[i]
                if i < int(max_sink_depth.magnitude / (2 * resolution.magnitude)):
                    p_allow = p_allow[:i * 2]
                    p_load = p_load[:i * 2]
                    depth = depth[:i * 2]
                break

        if sink_depth == -1.:
            print('Soil bearing capacity insuficient, solution does not converge within :' + str(max_sink_depth))

        return [p_allow.to('Pa'), p_load.to('Pa'), depth.to('m'), np.array([sink_depth.magnitude]) * ureg['m'],
                np.array([load.magnitude]) * ureg['Pa']]

    def solve_torque(self, depth, load):
        self.bot.depth = depth.copy()
        torque_req = self.bot.torque(load=load, layers=self.world.Layers)
        return torque_req

    def determine_max_torque(self):
        pass

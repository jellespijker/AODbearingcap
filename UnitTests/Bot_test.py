import unittest
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

from AOD.Material import *
from AOD.Bot import *
from AOD.Unit import *
from AOD.Model import *


class TestCylinderMethods(unittest.TestCase):
    def test_A_s(self):
        screw = Screw()
        screw.cylinder.d = 0.6 * ureg['m']
        screw.cylinder.l = 2. * ureg['m']
        depth = np.linspace(start=0, stop=screw.cylinder.d.magnitude, num=3) * ureg['m']
        a_s = screw.cylinder.A_s(depth) / screw.cylinder.A_t
        frac = np.array([0, 0.5, 1])
        np.testing.assert_array_equal(frac, a_s)

    def test_buoyancy(self):
        # Define the fluid en solid layers
        air = Air()
        water = Water()
        silt = Silt()
        layers = {'Air': air, 'Fluid': water, 'Soil': silt}

        screw = Screw()
        d = 0.6 * ureg['m']
        l = 2. * ureg['m']
        v = pi / 4 * d ** 2 * l
        screw.cylinder.l = l
        screw.cylinder.d = d
        depth = np.linspace(start=0, stop=screw.cylinder.d.magnitude, num=3) * ureg['m']
        f0 = v * g * (water.rho - air.rho)
        fr = 1 / 2 * v * g * (water.rho - air.rho) + 1 / 2 * v * g * (silt.rho_ins - air.rho)
        fd = v * g * (silt.rho_ins - air.rho)
        force = np.array([f0.magnitude, fr.magnitude, fd.magnitude]) * ureg['N']
        fb = screw.buoyancy(layers=layers, depth=depth)
        np.testing.assert_array_equal(np.round(fb, 2), np.round(force, 2))

    def test_b(self):
        screw = Screw()
        screw.cylinder.d = 0.6 * ureg['m']
        screw.cylinder.l = 2. * ureg['m']
        depth = np.linspace(start=0, stop=screw.cylinder.d.magnitude, num=3) * ureg['m']
        [b_acc, d_acc] = screw.cylinder.B_acc(depth)

        b_acc_exp = np.array([0., 0.6, 0.6]) * ureg['m']
        d_acc_exp = np.array([0., 0.235619, 0.471239]) * ureg['m']

        np.testing.assert_array_equal(np.round(b_acc, 1), np.round(b_acc_exp, 1))
        np.testing.assert_array_equal(np.round(d_acc, 4), np.round(d_acc_exp, 4))


class TestBotMethodes(unittest.TestCase):
    def test_dry_force(self):
        bot = Bot()
        bot.weight_dry = 5.1e3 * ureg['kg']
        force = 50013.915 * ureg['N']
        self.assertEqual(force, round(bot.F(), 3))

    def test_buoyancy_force(self):
        # Define the fluid en solid layers
        air = Air()
        water = Water()
        silt = Silt()
        layers = {'Air': air, 'Fluid': water, 'Soil': silt}

        bot = Bot()
        bot.weight_dry = 5.1e3 * ureg['kg']
        for s in bot.Screws:
            s.cylinder.l = 2. * ureg['m']
            s.cylinder.d = 0.6 * ureg['m']

        V = bot.Screw.cylinder.V

        f0 = bot.weight_dry * g - 2 * V * (water.rho - air.rho) * g
        fr = bot.weight_dry * g - (V * (water.rho - air.rho) * g + V * (silt.rho_ins - air.rho) * g)
        fd = bot.weight_dry * g - 2 * V * (silt.rho_ins - air.rho) * g
        fw = np.array([f0.magnitude, fr.magnitude, fd.magnitude]) * ureg['N']
        depth = np.linspace(start=0, stop=bot.Screw.cylinder.d.magnitude, num=3) * ureg['m']
        bot.depth = depth
        force = bot.F(forcetype=ForceType.Wet, layers=layers)
        np.testing.assert_array_equal(np.round(fw, 2), np.round(force, 2))

class TestModelMethods(unittest.TestCase):
    def test_solve(self):
        model = Model()
        #model.bot.weight_dry = 3.5e3 * ureg['kg']
        model.world.Layers['Soil'] = Packed_clay()
        depth = np.arange(start=0, stop=model.bot.Screw.cylinder.d.magnitude * 10, step=1.e-3) * ureg['m']
        [p_allow, p_load, sink_depth, load] = model.solve(depth=depth)
        print(sink_depth)
        print(load)
        plt.figure()
        plt.grid(True)
        plt.plot(sink_depth, load, '*')
        plt.plot(depth, p_allow)
        plt.plot(depth, p_load)
        plt.show()


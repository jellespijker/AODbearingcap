import unittest
import matplotlib.pyplot as plt
from AOD.Model import *
from AOD.Unit import *

class TestModelMethods(unittest.TestCase):
    def test_solve(self):
        model = Model()
        model.world.Layers['Soil'] = Sand()

        [p_allow, p_load, depth, sink_depth, load] = model.solve_sinkdepth()
        print(sink_depth)
        print(load)
        plt.figure()
        plt.grid(True)
        plt.plot(sink_depth, load, '*')
        plt.plot(depth, p_allow)
        plt.plot(depth, p_load)
        plt.show()
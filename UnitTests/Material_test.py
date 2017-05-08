import unittest
import matplotlib.pyplot as plt
from AOD.Material import *
from AOD.Unit import *

class TestSoilMethods(unittest.TestCase):
    def test_alpha(self):
        silt = Silt()
        self.assertEqual(silt.alpha(weight=0 * ureg['kg'], area=300 * ureg['m**2']), 0.95)
        self.assertEqual(silt.alpha(weight=2000 * ureg['kg'], area=1 * ureg['m**2']), 0.39)
        self.assertEqual(silt.alpha(weight=1. * ureg['kg'], area=2. * ureg['m**2']), 0.65)
        self.assertEqual(silt.alpha(weight=1. * ureg['kg'], area=2.2 * ureg['m**2']), 0.6727272727272727)



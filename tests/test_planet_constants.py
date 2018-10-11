########### Standard ###########
import unittest
import logging

########### Local ###########=
from orpytal import planet_constants

########### External ###########
import numpy as np

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)

class TestOrbit(unittest.TestCase):
    def test_json(self):
        moon = planet_constants.moon
        moon2 = planet_constants.CentralBody.from_dict(moon.to_dict())

        assert (moon.a == moon2.a)
        assert (moon.radius == moon2.radius)
        assert (moon.parent.a == moon2.parent.a)
        assert (moon.parent.period == moon2.parent.period)
        assert (moon.parent.parent.radius == moon2.parent.parent.radius)

if __name__ == '__main__':
    unittest.main()

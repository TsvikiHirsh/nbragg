import unittest
import numpy as np
from models import TransmissionModel
from lmfit import Parameters

# A mock cross-section function for testing
class MockCrossSection:
    def __init__(self, isotopes):
        self.isotopes = isotopes
    
    def __call__(self, E):
        # Return a dummy cross-section value that scales with energy
        return sum(self.isotopes.values()) * np.sqrt(E)

class TestTransmissionModel(unittest.TestCase):

    def setUp(self):
        """Set up the model with a mock cross section for testing."""
        isotopes = {"H-1": 0.5, "O-16": 0.3}
        self.cross_section = MockCrossSection(isotopes)
        self.model = TransmissionModel(self.cross_section)

    def test_initial_parameters(self):
        """Test the initial parameters are correctly set."""
        params = self.model.params
        self.assertAlmostEqual(params['H1'].value, 0.5)
        self.assertAlmostEqual(params['O16'].value, 0.3)
        self.assertFalse(params['H1'].vary)
        self.assertFalse(params['O16'].vary)
        self.assertFalse(params['n'].vary)

    def test_transmission_calculation(self):
        """Test the transmission function calculation."""
        E = np.array([1.0, 10.0, 100.0])
        T = self.model.transmission(E, thickness=1, n=0.01)
        self.assertEqual(T.shape, E.shape)
        self.assertTrue(np.all(T > 0))

    def test_background_varying(self):
        """Test the model with varying background parameters."""
        model_vary_bg = TransmissionModel(self.cross_section, vary_background=True)
        params = model_vary_bg.params
        self.assertTrue(params['b0'].vary)
        self.assertTrue(params['b1'].vary)
        self.assertTrue(params['b2'].vary)

if __name__ == '__main__':
    unittest.main()

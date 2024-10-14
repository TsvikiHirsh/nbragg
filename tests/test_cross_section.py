import unittest
import numpy as np
from cross_section import CrossSection, grab_from_endf

class TestCrossSection(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Set up a CrossSection object for Al-27 before any tests run."""
        cls.cross_section_al27 = CrossSection(isotopes={"Al-27": 1.0})

    def test_cross_section_al27(self):
        """Test CrossSection values for Al-27 at specific energies."""
        # Test at 1.0 MeV
        energy_1 = np.array([1.0])
        expected_xs_1 = self.cross_section_al27(energy_1)[0]
        self.assertAlmostEqual(expected_xs_1, 2.36756, delta=1e-5, 
                               msg=f"Expected cross section at 1.0 MeV: ~2.36756, but got {expected_xs_1}")

        # Test at 0.1 MeV
        energy_2 = np.array([0.1])
        expected_xs_2 = self.cross_section_al27(energy_2)[0]
        self.assertAlmostEqual(expected_xs_2, 5.30213310, delta=1e-4,
                               msg=f"Expected cross section at 0.1 MeV: ~5.30213, but got {expected_xs_2}")

        # Test at 10 MeV
        energy_3 = np.array([10.0])
        expected_xs_3 = self.cross_section_al27(energy_3)[0]
        self.assertAlmostEqual(expected_xs_3, 1.723, delta=1e-3,
                               msg=f"Expected cross section at 10 MeV: ~1.723, but got {expected_xs_3}")

    def test_cross_section_combination(self):
        """Test CrossSection values for a combination of isotopes."""
        isotopes = {"Al-27": 0.7, "O-16": 0.3}
        combined_cs = CrossSection(isotopes=isotopes)
        
        # Test at 1.0 MeV
        energy_1 = np.array([1.0])
        combined_xs_1 = combined_cs(energy_1)[0]
        self.assertAlmostEqual(combined_xs_1, 4.09538, delta=1e-3,
                               msg=f"Expected combined cross section at 1.0 MeV: ~4.09538, but got {combined_xs_1}")

        # Test at 0.1 MeV
        energy_2 = np.array([0.1])
        combined_xs_2 = combined_cs(energy_2)[0]
        self.assertAlmostEqual(combined_xs_2, 4.78696, delta=1e-3,
                               msg=f"Expected combined cross section at 0.1 MeV: ~4.78696, but got {combined_xs_2}")

        # Test at 10 MeV
        energy_3 = np.array([10.0])
        combined_xs_3 = combined_cs(energy_3)[0]
        self.assertAlmostEqual(combined_xs_3, 1.600516, delta=1e-3,
                               msg=f"Expected combined cross section at 10 MeV: ~1.600516, but got {combined_xs_3}")

    def test_cross_section_non_existing_isotope(self):
        """Test handling of non-existing isotopes."""
        with self.assertRaises(Exception):
            grab_from_endf("Xy-999")


if __name__ == "__main__":
    unittest.main()

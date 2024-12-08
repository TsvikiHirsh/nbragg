import unittest
import numpy as np
import nbragg

class TestMTEXToNCrystalConversion(unittest.TestCase):
    def setUp(self):
        # Path to your test CSV file
        self.csv_file = "tests/simple_components.csv"
        self.base_material = nbragg.materials["Fe_sg225_Iron-gamma.ncmat"]


    def test_first_phase_orientation(self):
        # Create CrossSection from MTEX data
        cs = nbragg.CrossSection().from_mtex(self.csv_file, self.base_material, short_name="γ")
        
        # Check the first phase (γ1)
        first_phase = cs.materials['γ1']
        
        # Normalized dir1 should be [0, 0, 1] (beam direction)
        expected_dir1 = [0, 0, 1.0]
        np.testing.assert_almost_equal(first_phase['dir1'], expected_dir1, decimal=7)
        
        # Normalized dir2 should be towards y-axis 
        expected_dir2 = [0, 1.0, 0]
        np.testing.assert_almost_equal(first_phase['dir2'], expected_dir2, decimal=7)
        
        # Check other properties
        self.assertEqual(first_phase['temp'], 300.0)
        self.assertEqual(first_phase['mos'], 10.0)
        self.assertAlmostEqual(first_phase['weight'], 1/7, places=7)

    def test_phases_object_creation(self):
        # Create CrossSection from MTEX data
        cs = nbragg.CrossSection().from_mtex(self.csv_file, self.base_material, short_name="γ")
        
        # Check phases object creation
        phases = cs.phases
        
        # Check number of phases
        self.assertEqual(len(phases), 7)
        
        # Check first phase details
        first_phase = phases['γ1']
        
        # Verify phase string format
        expected_prefix = 'Fe_sg225_Iron-gamma.ncmat;temp=300K;mos=10.0deg;dirtol=1.0deg;'
        self.assertTrue(first_phase.startswith(expected_prefix))
        
        # Check dir1 and dir2 parts
        self.assertIn('dir1=@crys_hkl:0.00000000,0.00000000,1.00000000@lab:0,0,1', first_phase)
        self.assertIn('dir2=@crys_hkl:0.00000000,1.00000000,0.00000000@lab:0,1,0', first_phase)

if __name__ == '__main__':
    unittest.main()
import unittest
import numpy as np
import nbragg

import unittest
import numpy as np
import nbragg
from nbragg.utils import materials as materials_dict

class TestCrossSection(unittest.TestCase):
    def test_cross_section_init_with_materials_dict(self):
        """Test initialization with dictionary of materials."""
        xs = nbragg.CrossSection(
            gamma=materials_dict["Fe_sg225_Iron-gamma.ncmat"],
            alpha="Fe_sg229_Iron-alpha.ncmat"
        )
        
        self.assertEqual(len(xs.materials), 2)
        
        # Check gamma material
        self.assertEqual(xs.materials['gamma']['mat'], 'Fe_sg225_Iron-gamma.ncmat')
        self.assertEqual(xs.materials['gamma']['temp'], 300.0)
        self.assertAlmostEqual(xs.materials['gamma']['weight'], 0.5)
        
        # Check alpha material
        self.assertEqual(xs.materials['alpha']['mat'], 'Fe_sg229_Iron-alpha.ncmat')
        self.assertEqual(xs.materials['alpha']['temp'], 300.0)
        self.assertAlmostEqual(xs.materials['alpha']['weight'], 0.5)

    def test_cross_section_init_with_custom_weights(self):
        """Test initialization with custom weights."""
        xs = nbragg.CrossSection({
            'gamma': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat', 
                'weight': 0.7
            },
            'alpha': {
                'mat': 'Fe_sg229_Iron-alpha.ncmat', 
                'weight': 0.3
            }
        })
        
        self.assertEqual(len(xs.materials), 2)
        
        self.assertEqual(xs.materials['gamma']['mat'], 'Fe_sg225_Iron-gamma.ncmat')
        self.assertAlmostEqual(xs.materials['gamma']['weight'], 0.7)
        
        self.assertEqual(xs.materials['alpha']['mat'], 'Fe_sg229_Iron-alpha.ncmat')
        self.assertAlmostEqual(xs.materials['alpha']['weight'], 0.3)

    def test_cross_section_init_with_total_weight(self):
        """Test initialization with total weight scaling."""
        xs = nbragg.CrossSection(
            gamma=materials_dict["Fe_sg225_Iron-gamma.ncmat"],
            alpha="Fe_sg229_Iron-alpha.ncmat", 
            total_weight=2.0
        )
        
        # Verify that individual material weights are scaled
        self.assertAlmostEqual(xs.materials['gamma']['weight'], 1.0)
        self.assertAlmostEqual(xs.materials['alpha']['weight'], 1.0)

    def test_cross_section_add_operator(self):
        """Test addition of two CrossSection objects."""
        xs1 = nbragg.CrossSection(gamma=materials_dict["Fe_sg225_Iron-gamma.ncmat"])
        xs2 = nbragg.CrossSection(alpha="Fe_sg229_Iron-alpha.ncmat")
        
        xs_combined = xs1 + xs2
        
        self.assertEqual(len(xs_combined.materials), 2)
        self.assertIn('gamma', xs_combined.materials)
        self.assertIn('alpha', xs_combined.materials)

    def test_cross_section_multiply_operator(self):
        """Test multiplication of a CrossSection by a scalar."""
        xs1 = nbragg.CrossSection(
            gamma=materials_dict["Fe_sg225_Iron-gamma.ncmat"], 
            total_weight=1.0
        )
        
        xs_scaled = xs1 * 2.0
        
        self.assertEqual(xs_scaled.total_weight, 2.0)

    def test_cross_section_with_orientation(self):
        """Test initialization with material orientation parameters."""
        xs = nbragg.CrossSection({
            'gamma': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'mos': 0.5,
                'dir1': [0, 0, 1],
                'dir2': [1, 0, 0],
                'theta': 45,
                'phi': 30
            }
        })
        
        self.assertEqual(xs.materials['gamma']['mos'], 0.5)
        self.assertEqual(xs.materials['gamma']['dir1'], [0, 0, 1])
        self.assertEqual(xs.materials['gamma']['dir2'], [1, 0, 0])
        self.assertEqual(xs.materials['gamma']['theta'], 45)
        self.assertEqual(xs.materials['gamma']['phi'], 30)

    def test_cross_section_with_temperature(self):
        """Test initialization with custom temperature."""
        xs = nbragg.CrossSection({
            'gamma': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'temp': 500
            }
        })
        
        self.assertEqual(xs.materials['gamma']['temp'], 500)

    def test_cross_section_nested_weight_normalization(self):
        """Test weight normalization with multiple materials."""
        xs = nbragg.CrossSection({
            'gamma1': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'weight': 0.3
            },
            'gamma2': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'weight': 0.7
            }
        })
        
        # Weights should sum to 1
        total_weight = sum(mat['weight'] for mat in xs.materials.values())
        self.assertAlmostEqual(total_weight, 1.0)

    def test_cross_section_with_string_material_ref(self):
        """Test initialization using string references to materials."""
        xs = nbragg.CrossSection(
            iron_gamma = 'Fe_sg225_Iron-gamma.ncmat',
            iron_alpha = 'Fe_sg229_Iron-alpha.ncmat'
        )
        
        self.assertEqual(xs.materials['iron_gamma']['mat'], 'Fe_sg225_Iron-gamma.ncmat')
        self.assertEqual(xs.materials['iron_alpha']['mat'], 'Fe_sg229_Iron-alpha.ncmat')


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
        expected_prefix = 'Fe_sg225_Iron-gamma.nbragg;temp=300K;mos=10.0deg;dirtol=1.0deg;'
        self.assertTrue(first_phase.startswith(expected_prefix))
        
        # Check dir1 and dir2 parts
        self.assertIn('dir1=@crys_hkl:0.00000000,0.00000000,1.00000000@lab:0,0,1', first_phase)
        self.assertIn('dir2=@crys_hkl:0.00000000,1.00000000,0.00000000@lab:0,1,0', first_phase)

if __name__ == '__main__':
    unittest.main()
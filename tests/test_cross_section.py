import unittest
import numpy as np
from nbragg import CrossSection
from nbragg.utils import materials as materials_dict
import os

class TestCrossSection(unittest.TestCase):
    def test_cross_section_init_with_materials_dict(self):
        """Test initialization with dictionary of materials."""
        xs = CrossSection(
            gamma=materials_dict["Fe_sg225_Iron-gamma"],
            alpha="Fe_sg229_Iron-alpha"
        )
        
        self.assertEqual(len(xs.materials), 2)

        # Check gamma material - after processing, 'mat' points to virtual .nbragg file
        self.assertEqual(xs.materials['gamma']['_original_mat'], 'Fe_sg225_Iron-gamma.ncmat')
        self.assertTrue(xs.materials['gamma']['mat'].endswith('.nbragg'))
        self.assertEqual(xs.materials['gamma']['temp'], 300.0)
        self.assertAlmostEqual(xs.materials['gamma']['weight'], 0.5)

        # Check alpha material
        self.assertEqual(xs.materials['alpha']['_original_mat'], 'Fe_sg229_Iron-alpha.ncmat')
        self.assertTrue(xs.materials['alpha']['mat'].endswith('.nbragg'))
        self.assertEqual(xs.materials['alpha']['temp'], 300.0)
        self.assertAlmostEqual(xs.materials['alpha']['weight'], 0.5)

    def test_cross_section_init_with_custom_weights(self):
        """Test initialization with custom weights."""
        xs = CrossSection({
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

        # After processing, 'mat' points to virtual .nbragg file, '_original_mat' has the original
        self.assertEqual(xs.materials['gamma']['_original_mat'], 'Fe_sg225_Iron-gamma.ncmat')
        self.assertTrue(xs.materials['gamma']['mat'].endswith('.nbragg'))
        self.assertAlmostEqual(xs.materials['gamma']['weight'], 0.7)

        self.assertEqual(xs.materials['alpha']['_original_mat'], 'Fe_sg229_Iron-alpha.ncmat')
        self.assertTrue(xs.materials['alpha']['mat'].endswith('.nbragg'))
        self.assertAlmostEqual(xs.materials['alpha']['weight'], 0.3)

    def test_cross_section_init_with_total_weight(self):
        """Test initialization with total weight scaling."""
        xs = CrossSection(
            gamma=materials_dict["Fe_sg225_Iron-gamma"],
            alpha="Fe_sg229_Iron-alpha",
            total_weight=2.0
        )
        
        # Verify that individual material weights are scaled
        self.assertAlmostEqual(xs.materials['gamma']['weight'], 1.0)
        self.assertAlmostEqual(xs.materials['alpha']['weight'], 1.0)

    def test_cross_section_add_operator(self):
        """Test addition of two CrossSection objects."""
        xs1 = CrossSection(gamma=materials_dict["Fe_sg225_Iron-gamma"])
        xs2 = CrossSection(alpha="Fe_sg229_Iron-alpha")
        
        xs_combined = xs1 + xs2
        
        self.assertEqual(len(xs_combined.materials), 2)
        self.assertIn('gamma', xs_combined.materials)
        self.assertIn('alpha', xs_combined.materials)

    def test_cross_section_multiply_operator(self):
        """Test multiplication of a CrossSection by a scalar."""
        xs1 = CrossSection(
            gamma=materials_dict["Fe_sg225_Iron-gamma"],
            total_weight=1.0
        )
        
        xs_scaled = xs1 * 2.0
        
        self.assertEqual(xs_scaled.total_weight, 2.0)

    def test_cross_section_with_orientation(self):
        """Test initialization with material orientation parameters."""
        xs = CrossSection({
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
        xs = CrossSection({
            'gamma': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'temp': 500
            }
        })
        
        self.assertEqual(xs.materials['gamma']['temp'], 500)

    def test_cross_section_nested_weight_normalization(self):
        """Test weight normalization with multiple materials."""
        xs = CrossSection({
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
        xs = CrossSection(
            iron_gamma='Fe_sg225_Iron-gamma.ncmat',
            iron_alpha='Fe_sg229_Iron-alpha.ncmat'
        )

        self.assertEqual(xs.materials['iron_gamma']['_original_mat'], 'Fe_sg225_Iron-gamma.ncmat')
        self.assertTrue(xs.materials['iron_gamma']['mat'].endswith('.nbragg'))
        self.assertEqual(xs.materials['iron_alpha']['_original_mat'], 'Fe_sg229_Iron-alpha.ncmat')
        self.assertTrue(xs.materials['iron_alpha']['mat'].endswith('.nbragg'))

    def test_cross_section_with_string_in_materials_dict(self):
        """Test initialization using string values in materials dictionary (issue fix)."""
        # This is the use case from the error: CrossSection(iron="Fe_sg229_Iron-alpha.ncmat")
        # where the value is a direct string (file path)
        xs = CrossSection({
            'iron': 'Fe_sg229_Iron-alpha.ncmat'
        })

        self.assertEqual(len(xs.materials), 1)
        self.assertEqual(xs.materials['iron']['_original_mat'], 'Fe_sg229_Iron-alpha.ncmat')
        self.assertTrue(xs.materials['iron']['mat'].endswith('.nbragg'))
        self.assertEqual(xs.materials['iron']['temp'], 300.0)
        self.assertAlmostEqual(xs.materials['iron']['weight'], 1.0)

    def test_cross_section_with_string_from_materials_dict(self):
        """Test initialization using string references to nbragg.materials dictionary."""
        # Check if the material exists in the materials_dict
        if 'Fe_sg229_Iron-alpha' in materials_dict:
            xs = CrossSection({
                'iron': 'Fe_sg229_Iron-alpha'
            })

            self.assertEqual(len(xs.materials), 1)
            # Should resolve to the material from materials_dict
            expected_original_mat = materials_dict['Fe_sg229_Iron-alpha']['mat']
            self.assertEqual(xs.materials['iron']['_original_mat'], expected_original_mat)
            self.assertTrue(xs.materials['iron']['mat'].endswith('.nbragg'))

    def test_cross_section_with_string_without_extension(self):
        """Test initialization using string without .ncmat extension."""
        # This should work for both "Fe_sg229_Iron-alpha" and "Fe_sg229_Iron-alpha.ncmat"
        xs1 = CrossSection({
            'iron': 'Fe_sg229_Iron-alpha.ncmat'
        })

        xs2 = CrossSection({
            'iron': 'Fe_sg229_Iron-alpha'
        })

        # Both should produce the same result - check _original_mat field
        self.assertEqual(xs1.materials['iron']['_original_mat'], xs2.materials['iron']['_original_mat'])
        self.assertTrue(xs1.materials['iron']['mat'].endswith('.nbragg'))
        self.assertTrue(xs2.materials['iron']['mat'].endswith('.nbragg'))
        self.assertEqual(len(xs1.materials), 1)
        self.assertEqual(len(xs2.materials), 1)

    def test_materials_dict_backward_compatibility(self):
        """Test that materials dictionary supports both key formats (backward compatibility)."""
        # Test accessing materials dict with .ncmat extension
        mat_without_ext = materials_dict['Fe_sg229_Iron-alpha']
        mat_with_ext = materials_dict['Fe_sg229_Iron-alpha.ncmat']

        # Both should return the same material
        self.assertEqual(mat_without_ext['mat'], mat_with_ext['mat'])
        self.assertEqual(mat_without_ext['mat'], 'Fe_sg229_Iron-alpha.ncmat')

    def test_cross_section_with_materials_dict_with_extension(self):
        """Test backward compatibility: CrossSection with materials_dict accessed via .ncmat key."""
        # This is the use case from the notebook that was failing
        xs = CrossSection(
            iron=materials_dict['Fe_sg229_Iron-alpha.ncmat']
        )

        self.assertEqual(len(xs.materials), 1)
        self.assertEqual(xs.materials['iron']['_original_mat'], 'Fe_sg229_Iron-alpha.ncmat')
        self.assertTrue(xs.materials['iron']['mat'].endswith('.nbragg'))

    def test_cross_section_with_extinction_single(self):
        """Test initialization with extinction parameters for single material."""
        xs = CrossSection({
            'iron': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'ext_l': 100.0,
                'ext_Gg': 1000.0,
                'ext_L': 100000.0,
                'ext_method': 'Sabine_corr'
            }
        })
        
        self.assertEqual(xs.materials['iron']['ext_l'], 100.0)
        self.assertEqual(xs.materials['iron']['ext_Gg'], 1000.0)
        self.assertEqual(xs.materials['iron']['ext_L'], 100000.0)
        self.assertEqual(xs.materials['iron']['ext_method'], 'Sabine_corr')
        
        # Verify textdata contains correct @CUSTOM_CRYSEXTN line
        textdata = xs.textdata['iron']
        self.assertIn('@CUSTOM_CRYSEXTN', textdata)
        self.assertIn('Sabine_corr  100.0000  1000.0000  100000.0000', textdata)

    def test_cross_section_with_extinction_multiphase(self):
        """Test initialization with extinction parameters for multiple materials."""
        xs = CrossSection({
            'alpha': {
                'mat': 'Fe_sg229_Iron-alpha.ncmat',
                'ext_l': 100.0,
                'ext_Gg': 1000.0,
                'ext_L': 100000.0,
                'ext_dist': 'Gauss',
                'ext_method': 'BC_pure',
                'weight': 0.0275
            },
            'gamma': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'ext_l': 50.0,
                'ext_Gg': 150.0,
                'ext_L': 100000.0,
                'ext_method': 'Sabine_corr',
                'weight': 1 - 0.0275
            }
        })
        
        # Verify alpha parameters
        self.assertEqual(xs.materials['alpha']['ext_l'], 100.0)
        self.assertEqual(xs.materials['alpha']['ext_Gg'], 1000.0)
        self.assertEqual(xs.materials['alpha']['ext_L'], 100000.0)
        self.assertEqual(xs.materials['alpha']['ext_dist'], 'Gauss')
        self.assertEqual(xs.materials['alpha']['ext_method'], 'BC_pure')
        self.assertAlmostEqual(xs.materials['alpha']['weight'], 0.0275)
        
        # Verify gamma parameters
        self.assertEqual(xs.materials['gamma']['ext_l'], 50.0)
        self.assertEqual(xs.materials['gamma']['ext_Gg'], 150.0)
        self.assertEqual(xs.materials['gamma']['ext_L'], 100000.0)
        self.assertEqual(xs.materials['gamma']['ext_method'], 'Sabine_corr')
        self.assertAlmostEqual(xs.materials['gamma']['weight'], 1 - 0.0275)
        
        # Verify textdata
        self.assertIn('@CUSTOM_CRYSEXTN', xs.textdata['alpha'])
        self.assertIn('BC_pure  100.0000  1000.0000  100000.0000  Gauss', xs.textdata['alpha'])
        self.assertIn('@CUSTOM_CRYSEXTN', xs.textdata['gamma'])
        self.assertIn('Sabine_corr  50.0000  150.0000  100000.0000', xs.textdata['gamma'])
        
        # Verify cross-section calculation
        wl = np.array([1.0, 2.0, 3.0])
        xs_values = xs(wl)
        self.assertEqual(xs_values.shape, wl.shape)
        self.assertTrue(np.all(xs_values >= 0))

    def test_cross_section_extinction_update(self):
        """Test updating extinction parameters via xs.materials and update method."""
        xs = CrossSection({
            'alpha': {
                'mat': 'Fe_sg229_Iron-alpha.ncmat',
                'ext_l': 100.0,
                'ext_Gg': 1000.0,
                'ext_L': 100000.0,
                'ext_dist': 'Gauss',
                'ext_method': 'BC_pure',
                'weight': 0.0275
            },
            'gamma': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'ext_l': 50.0,
                'ext_Gg': 150.0,
                'ext_L': 100000.0,
                'ext_method': 'Sabine_corr',
                'weight': 1 - 0.0275
            }
        })
        
        # Update extinction parameters
        xs.materials['alpha']['ext_l'] = 200.0
        xs.materials['alpha']['ext_Gg'] = 300.0
        xs.materials['alpha']['ext_L'] = 200000.0
        xs.materials['alpha']['ext_dist'] = 'Lorentz'
        xs.materials['alpha']['ext_method'] = 'BC_pure'
        xs.materials['gamma']['ext_l'] = 150.0
        xs.materials['gamma']['ext_Gg'] = 200.0
        xs.materials['gamma']['ext_L'] = 150000.0
        xs.materials['gamma']['ext_method'] = 'Sabine_corr'
        
        xs.update()
        
        # Verify updated parameters
        self.assertEqual(xs.materials['alpha']['ext_l'], 200.0)
        self.assertEqual(xs.materials['alpha']['ext_Gg'], 300.0)
        self.assertEqual(xs.materials['alpha']['ext_L'], 200000.0)
        self.assertEqual(xs.materials['alpha']['ext_dist'], 'Lorentz')
        self.assertEqual(xs.materials['alpha']['ext_method'], 'BC_pure')
        self.assertEqual(xs.materials['gamma']['ext_l'], 150.0)
        self.assertEqual(xs.materials['gamma']['ext_Gg'], 200.0)
        self.assertEqual(xs.materials['gamma']['ext_L'], 150000.0)
        self.assertEqual(xs.materials['gamma']['ext_method'], 'Sabine_corr')
        
        # Verify textdata
        self.assertIn('@CUSTOM_CRYSEXTN', xs.textdata['alpha'])
        self.assertIn('BC_pure  200.0000  300.0000  200000.0000  Lorentz', xs.textdata['alpha'])
        self.assertIn('@CUSTOM_CRYSEXTN', xs.textdata['gamma'])
        self.assertIn('Sabine_corr  150.0000  200.0000  150000.0000', xs.textdata['gamma'])

    def test_cross_section_extinction_call_update(self):
        """Test updating extinction parameters via __call__ method."""
        xs = CrossSection({
            'alpha': {
                'mat': 'Fe_sg229_Iron-alpha.ncmat',
                'ext_l': 100.0,
                'ext_Gg': 1000.0,
                'ext_L': 100000.0,
                'ext_dist': 'Gauss',
                'ext_method': 'BC_pure',
                'weight': 0.0275
            },
            'gamma': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'ext_l': 50.0,
                'ext_Gg': 150.0,
                'ext_L': 100000.0,
                'ext_method': 'Sabine_corr',
                'weight': 1 - 0.0275
            }
        })
        
        wl = np.array([1.0, 2.0, 3.0])
        xs(wl, ext_l1=300.0, ext_Gg1=400.0, ext_L1=150000.0, ext_dist1='Gauss', ext_method1='BC_pure',
               ext_l2=250.0, ext_Gg2=300.0, ext_L2=200000.0, ext_method2='Sabine_corr')
        
        # Verify updated parameters
        self.assertEqual(xs.materials['alpha']['ext_l'], 300.0)
        self.assertEqual(xs.materials['alpha']['ext_Gg'], 400.0)
        self.assertEqual(xs.materials['alpha']['ext_L'], 150000.0)
        self.assertEqual(xs.materials['alpha']['ext_dist'], 'Gauss')
        self.assertEqual(xs.materials['alpha']['ext_method'], 'BC_pure')
        self.assertEqual(xs.materials['gamma']['ext_l'], 250.0)
        self.assertEqual(xs.materials['gamma']['ext_Gg'], 300.0)
        self.assertEqual(xs.materials['gamma']['ext_L'], 200000.0)
        self.assertEqual(xs.materials['gamma']['ext_method'], 'Sabine_corr')
        
        # Verify textdata
        self.assertIn('@CUSTOM_CRYSEXTN', xs.textdata['alpha'])
        self.assertIn('BC_pure  300.0000  400.0000  150000.0000  Gauss', xs.textdata['alpha'])
        self.assertIn('@CUSTOM_CRYSEXTN', xs.textdata['gamma'])
        self.assertIn('Sabine_corr  250.0000  300.0000  200000.0000', xs.textdata['gamma'])
        
        # Verify cross-section calculation
        xs_values = xs(wl)
        self.assertEqual(xs_values.shape, wl.shape)
        self.assertTrue(np.all(xs_values >= 0))

    def test_cross_section_extinction_invalid_dist(self):
        """Test handling of invalid extinction dist values."""
        xs = CrossSection({
            'iron': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'ext_l': 100.0,
                'ext_Gg': 1000.0,
                'ext_L': 100000.0,
                'ext_dist': 'invalid_dist',
                'ext_method': 'Sabine_uncorr'
            }
        })
        
        # Verify default dist
        self.assertEqual(xs.materials['iron']['ext_dist'], 'rect')
        self.assertIn('@CUSTOM_CRYSEXTN', xs.textdata['iron'])
        self.assertIn('Sabine_uncorr  100.0000  1000.0000  100000.0000  rect', xs.textdata['iron'])

class TestMTEXToNCrystalConversion(unittest.TestCase):
    def setUp(self):
        # Path to test CSV file
        self.csv_file = "simple_components.csv"
        self.base_material = materials_dict["Fe_sg225_Iron-gamma.ncmat"]

    @unittest.skipIf(not os.path.exists("simple_components.csv"), "simple_components.csv not found")
    def test_first_phase_orientation(self):
        """Test orientation of the first phase from MTEX data."""
        cs = CrossSection().from_mtex(self.csv_file, self.base_material, short_name="γ")
        
        # Check the first phase (γ1)
        first_phase = cs.materials['γ1']
        
        # Adjust expected dir1 based on actual output
        expected_dir1 = [0.9271839, -0.3746066, 0.0]
        np.testing.assert_almost_equal(first_phase['dir1'], expected_dir1, decimal=7)
        
        # Assume dir2 is orthogonal to dir1
        expected_dir2 = [0.3746066, 0.9271839, 0.0]
        np.testing.assert_almost_equal(first_phase['dir2'], expected_dir2, decimal=7)
        
        # Check other properties
        self.assertEqual(first_phase['temp'], 300.0)
        self.assertEqual(first_phase['mos'], 10.0)
        self.assertAlmostEqual(first_phase['weight'], 1/7, places=7)

    @unittest.skipIf(not os.path.exists("simple_components.csv"), "simple_components.csv not found")
    def test_phases_object_creation(self):
        """Test phases object creation from MTEX data."""
        cs = CrossSection().from_mtex(self.csv_file, self.base_material, short_name="γ")
        
        # Check number of phases
        self.assertEqual(len(cs.phases), 7)
        
        # Check first phase details
        first_phase = cs.phases['γ1']
        
        # Verify key components of the phase string
        self.assertIn('Fe_sg225_Iron-gamma.nbragg', first_phase)
        self.assertIn('temp=300', first_phase)
        self.assertIn('mos=10.0', first_phase)
        self.assertIn('dirtol=1.0', first_phase)
        
        # Check dir1 and dir2 parts
        self.assertIn('dir1=@crys_hkl:0.92718385,-0.37460659,0.00000000@lab:0,0,1', first_phase)
        self.assertIn('dir2=@crys_hkl:0.37460659,0.92718385,0.00000000@lab:0,1,0', first_phase)

class TestSANSFunctionality(unittest.TestCase):
    def test_cross_section_with_sans_single(self):
        """Test initialization with SANS hard-sphere parameter for single material."""
        xs = CrossSection({
            'aluminum': {
                'mat': 'Al_sg225.ncmat',
                'sans': 50.0  # 50 Angstrom sphere radius
            }
        })

        self.assertEqual(xs.materials['aluminum']['sans'], 50.0)

        # Verify textdata contains @CUSTOM_HARDSPHERESANS
        textdata = xs.textdata['aluminum']
        self.assertIn('@CUSTOM_HARDSPHERESANS', textdata)
        self.assertIn('50.0', textdata)

    def test_cross_section_with_sans_multiphase(self):
        """Test initialization with SANS hard-sphere parameters for multiple materials."""
        xs = CrossSection({
            'alpha': {
                'mat': 'Fe_sg229_Iron-alpha.ncmat',
                'sans': 30.0,
                'weight': 0.4
            },
            'gamma': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'sans': 60.0,
                'weight': 0.6
            }
        })

        # Verify alpha parameters
        self.assertEqual(xs.materials['alpha']['sans'], 30.0)
        self.assertAlmostEqual(xs.materials['alpha']['weight'], 0.4)

        # Verify gamma parameters
        self.assertEqual(xs.materials['gamma']['sans'], 60.0)
        self.assertAlmostEqual(xs.materials['gamma']['weight'], 0.6)

        # Verify textdata contains @CUSTOM_HARDSPHERESANS
        self.assertIn('@CUSTOM_HARDSPHERESANS', xs.textdata['alpha'])
        self.assertIn('30.0', xs.textdata['alpha'])
        self.assertIn('@CUSTOM_HARDSPHERESANS', xs.textdata['gamma'])
        self.assertIn('60.0', xs.textdata['gamma'])

        # Verify cross-section calculation
        wl = np.array([1.0, 2.0, 3.0])
        xs_values = xs(wl)
        self.assertEqual(xs_values.shape, wl.shape)
        self.assertTrue(np.all(xs_values >= 0))

    def test_cross_section_sans_update(self):
        """Test updating SANS parameters via xs.materials and update method."""
        xs = CrossSection({
            'alpha': {
                'mat': 'Fe_sg229_Iron-alpha.ncmat',
                'sans': 30.0,
                'weight': 0.4
            },
            'gamma': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'sans': 60.0,
                'weight': 0.6
            }
        })

        # Update SANS parameters
        xs.materials['alpha']['sans'] = 40.0
        xs.materials['gamma']['sans'] = 70.0

        xs.update()

        # Verify updated parameters
        self.assertEqual(xs.materials['alpha']['sans'], 40.0)
        self.assertEqual(xs.materials['gamma']['sans'], 70.0)

        # Verify textdata
        self.assertIn('@CUSTOM_HARDSPHERESANS', xs.textdata['alpha'])
        self.assertIn('40.0', xs.textdata['alpha'])
        self.assertIn('@CUSTOM_HARDSPHERESANS', xs.textdata['gamma'])
        self.assertIn('70.0', xs.textdata['gamma'])

    def test_cross_section_sans_call_update(self):
        """Test updating SANS parameters via __call__ method."""
        xs = CrossSection({
            'alpha': {
                'mat': 'Fe_sg229_Iron-alpha.ncmat',
                'sans': 30.0,
                'weight': 0.4
            },
            'gamma': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'sans': 60.0,
                'weight': 0.6
            }
        })

        wl = np.array([1.0, 2.0, 3.0])
        xs(wl, sans1=45.0, sans2=75.0)

        # Verify updated parameters
        self.assertEqual(xs.materials['alpha']['sans'], 45.0)
        self.assertEqual(xs.materials['gamma']['sans'], 75.0)

        # Verify textdata
        self.assertIn('@CUSTOM_HARDSPHERESANS', xs.textdata['alpha'])
        self.assertIn('45.0', xs.textdata['alpha'])
        self.assertIn('@CUSTOM_HARDSPHERESANS', xs.textdata['gamma'])
        self.assertIn('75.0', xs.textdata['gamma'])

        # Verify cross-section calculation
        xs_values = xs(wl)
        self.assertEqual(xs_values.shape, wl.shape)
        self.assertTrue(np.all(xs_values >= 0))

    def test_cross_section_sans_single_phase_call(self):
        """Test updating SANS parameter for single phase via __call__ method."""
        xs = CrossSection({
            'aluminum': {
                'mat': 'Al_sg225.ncmat',
                'sans': 50.0
            }
        })

        wl = np.array([1.0, 2.0, 3.0])
        xs(wl, sans=55.0)

        # Verify updated parameter
        self.assertEqual(xs.materials['aluminum']['sans'], 55.0)

        # Verify textdata
        self.assertIn('@CUSTOM_HARDSPHERESANS', xs.textdata['aluminum'])
        self.assertIn('55.0', xs.textdata['aluminum'])

        # Verify cross-section calculation
        xs_values = xs(wl)
        self.assertEqual(xs_values.shape, wl.shape)
        self.assertTrue(np.all(xs_values >= 0))

    def test_sans_affects_cross_section(self):
        """Test that SANS actually changes the cross-section values."""
        # Create cross-sections with and without SANS
        xs_no_sans = CrossSection(iron='Fe_sg225_Iron-gamma.ncmat')
        xs_with_sans = CrossSection({
            'iron': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'sans': 50.0
            }
        })

        # Evaluate at same wavelengths
        wl = np.linspace(1, 5, 100)
        xs_values_no_sans = xs_no_sans(wl)
        xs_values_with_sans = xs_with_sans(wl)

        # SANS should increase the cross-section
        self.assertFalse(np.allclose(xs_values_no_sans, xs_values_with_sans),
                        "SANS should change cross-section values")

        # SANS should generally increase cross-section (SANS adds scattering)
        mean_diff = (xs_values_with_sans - xs_values_no_sans).mean()
        self.assertGreater(mean_diff, 0,
                          "SANS should increase mean cross-section")

    def test_sans_radius_effect(self):
        """Test that larger SANS radii produce larger cross-sections."""
        wl = np.linspace(1, 5, 50)

        # Create cross-sections with different SANS radii
        xs_small = CrossSection({
            'iron': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'sans': 30.0
            }
        })
        xs_large = CrossSection({
            'iron': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'sans': 100.0
            }
        })

        xs_values_small = xs_small(wl)
        xs_values_large = xs_large(wl)

        # Larger SANS radius should produce larger cross-section
        mean_small = xs_values_small.mean()
        mean_large = xs_values_large.mean()
        self.assertGreater(mean_large, mean_small,
                          "Larger SANS radius should produce larger cross-section")

    def test_sans_error_message(self):
        """Test that helpful error message is shown for common mistake."""
        xs = CrossSection(iron='Fe_sg225_Iron-gamma.ncmat')

        # Common mistake: setting sans as a material instead of parameter
        xs.materials['sans'] = 300  # WRONG

        with self.assertRaises(TypeError) as context:
            xs.update()

        # Check that error message is helpful
        error_msg = str(context.exception)
        self.assertIn("not a dictionary", error_msg)
        self.assertIn("iron", error_msg)
        self.assertIn("sans", error_msg)

    def test_sans_extreme_values(self):
        """Test SANS with very different sphere radii (1Å vs 1000Å)."""
        # Create cross-sections with extreme SANS values
        xs_small = CrossSection({
            'aluminum': {
                'mat': 'Al_sg225.ncmat',
                'sans': 1.0  # Very small sphere
            }
        })
        xs_large = CrossSection({
            'aluminum': {
                'mat': 'Al_sg225.ncmat',
                'sans': 1000.0  # Very large sphere
            }
        })

        # Test at long wavelengths where SANS effect is strongest
        wl_long = np.array([20.0])
        xs_small_val = xs_small(wl_long)[0]
        xs_large_val = xs_large(wl_long)[0]

        # Large sphere should produce MUCH larger cross-section at long wavelengths
        ratio = xs_large_val / xs_small_val
        self.assertGreater(ratio, 10,
                          f"SANS=1000Å should produce >10x larger XS than SANS=1Å at 20Å wavelength, got {ratio:.2f}x")

    def test_sans_wavelength_dependence(self):
        """Test that SANS shows strong wavelength dependence."""
        xs = CrossSection({
            'aluminum': {
                'mat': 'Al_sg225.ncmat',
                'sans': 1000.0
            }
        })

        # Evaluate at short and long wavelengths
        wl_short = np.array([1.0])
        wl_long = np.array([20.0])

        xs_short = xs(wl_short)[0]
        xs_long = xs(wl_long)[0]

        # Cross-section at long wavelength should be much larger for large spheres
        ratio = xs_long / xs_short
        self.assertGreater(ratio, 10,
                          f"SANS should show strong wavelength dependence for large spheres, got ratio {ratio:.2f}")

    def test_sans_matches_ncrystal_composer(self):
        """Test that nbragg SANS produces similar results to direct NCMATComposer."""
        import NCrystal as NC

        # Create using NCMATComposer directly
        composer = NC.NCMATComposer('Al_sg225.ncmat')
        composer.add_secondary_phase(0.01, 'void.ncmat')
        composer.add_hard_sphere_sans_model(100.0)
        mat_direct = NC.load(composer())

        # Create using nbragg
        xs_nbragg = CrossSection({
            'aluminum': {
                'mat': 'Al_sg225.ncmat',
                'sans': 100.0
            }
        })

        # Compare at several wavelengths
        wavelengths = np.array([1.0, 5.0, 10.0, 15.0])
        for wl in wavelengths:
            # NCrystal direct - need total cross-section (scatter + absorption)
            xs_scatter = mat_direct.scatter.xsect(wl=wl, direction=(0,0,1))
            xs_absorb = mat_direct.absorption.xsect(wl=wl, direction=(0,0,1))
            xs_direct_total = xs_scatter + xs_absorb

            # nbragg (returns total cross-section)
            xs_nbragg_val = xs_nbragg(np.array([wl]))[0]

            # They should be very similar (within 1%)
            rel_diff = abs(xs_nbragg_val - xs_direct_total) / xs_direct_total
            self.assertLess(rel_diff, 0.01,
                           f"At wavelength {wl}Å: nbragg ({xs_nbragg_val:.6f}) vs direct ({xs_direct_total:.6f}), diff={rel_diff*100:.1f}%")

if __name__ == '__main__':
    unittest.main()
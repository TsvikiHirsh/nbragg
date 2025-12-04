"""Test suite for CrossSection.from_mtex method"""
import unittest
import os
import tempfile
import pandas as pd
import nbragg


class TestFromMTEX(unittest.TestCase):
    """Test the CrossSection.from_mtex class method"""

    @classmethod
    def setUpClass(cls):
        """Create a temporary CSV file for testing"""
        cls.temp_dir = tempfile.mkdtemp()
        cls.csv_file = os.path.join(cls.temp_dir, 'test_orientations.csv')

        # Create test data matching the expected format
        test_data = {
            'alpha_mtex': [289.72674177, 347.55678233],
            'beta_mtex': [35.96254728, 45.14846774],
            'gamma_mtex': [93.53748033, 357.8898998],
            'volume_mtex': [0.4, 0.3],  # Total < 1 to test powder phase
            'xh': [2.64315179, 1.94886872],
            'xk': [-0.94753028, 0.68948715],
            'xl': [0.56783894, 1.98318538],
            'yh': [1.09976362, -0.53804307],
            'yk': [2.11879441, 2.77948231],
            'yl': [-1.58358338, -0.43760006],
            'fwhm': [12.587233, 17.831007]
        }

        df = pd.DataFrame(test_data)
        df.to_csv(cls.csv_file, index=False)

    @classmethod
    def tearDownClass(cls):
        """Clean up temporary files"""
        import shutil
        if os.path.exists(cls.temp_dir):
            shutil.rmtree(cls.temp_dir)

    def test_from_mtex_with_dict_and_short_name(self):
        """Test from_mtex with material dict and explicit short_name"""
        material = nbragg.materials["Fe_sg225_Iron-gamma.ncmat"]
        xs = nbragg.CrossSection.from_mtex(
            self.csv_file,
            material=material,
            short_name="gamma"
        )

        # Check that CrossSection was created
        self.assertIsInstance(xs, nbragg.CrossSection)

        # Check that we have oriented phases + powder phase
        self.assertIn('gamma0', xs.materials)
        self.assertIn('gamma1', xs.materials)
        self.assertIn('gamma_powder', xs.materials)

        # Check weights sum to 1
        total_weight = sum(mat['weight'] for mat in xs.materials.values())
        self.assertAlmostEqual(total_weight, 1.0, places=6)

        # Check that oriented phases have orientation parameters
        self.assertIsNotNone(xs.materials['gamma0']['dir1'])
        self.assertIsNotNone(xs.materials['gamma0']['dir2'])
        self.assertIsNotNone(xs.materials['gamma0']['mos'])

        # Check that powder phase has no orientation
        self.assertIsNone(xs.materials['gamma_powder']['dir1'])
        self.assertIsNone(xs.materials['gamma_powder']['dir2'])
        self.assertIsNone(xs.materials['gamma_powder']['mos'])

    def test_from_mtex_with_string_material(self):
        """Test from_mtex with material as string (filename)"""
        xs = nbragg.CrossSection.from_mtex(
            self.csv_file,
            material="Fe_sg225_Iron-gamma.ncmat",
            short_name="gamma"
        )

        # Check that CrossSection was created
        self.assertIsInstance(xs, nbragg.CrossSection)

        # Check that we have phases
        self.assertIn('gamma0', xs.materials)
        self.assertIn('gamma1', xs.materials)

        # Check that original mat field is preserved
        self.assertEqual(xs.materials['gamma0']['_original_mat'], "Fe_sg225_Iron-gamma.ncmat")

    def test_from_mtex_without_short_name(self):
        """Test from_mtex without short_name (should auto-generate)"""
        material = nbragg.materials["Fe_sg225_Iron-gamma.ncmat"]
        xs = nbragg.CrossSection.from_mtex(
            self.csv_file,
            material=material
        )

        # Check that CrossSection was created
        self.assertIsInstance(xs, nbragg.CrossSection)

        # Check that phases were created with auto-generated name
        # Should extract "Iron-gamma" from "Fe_sg225_Iron-gamma.ncmat"
        self.assertIn('Iron-gamma0', xs.materials)
        self.assertIn('Iron-gamma1', xs.materials)
        self.assertIn('Iron-gamma_powder', xs.materials)

    def test_from_mtex_without_short_name_string_material(self):
        """Test from_mtex with string material and no short_name"""
        xs = nbragg.CrossSection.from_mtex(
            self.csv_file,
            material="Fe_sg225_Iron-gamma.ncmat"
        )

        # Check that CrossSection was created
        self.assertIsInstance(xs, nbragg.CrossSection)

        # Should extract name from filename
        self.assertIn('Iron-gamma0', xs.materials)
        self.assertIn('Iron-gamma1', xs.materials)

    def test_from_mtex_no_powder_phase(self):
        """Test from_mtex with powder_phase=False"""
        material = nbragg.materials["Fe_sg225_Iron-gamma.ncmat"]
        xs = nbragg.CrossSection.from_mtex(
            self.csv_file,
            material=material,
            short_name="gamma",
            powder_phase=False
        )

        # Check that powder phase was NOT created
        self.assertNotIn('gamma_powder', xs.materials)

        # But oriented phases should still exist
        self.assertIn('gamma0', xs.materials)
        self.assertIn('gamma1', xs.materials)

    def test_from_mtex_weight_normalization(self):
        """Test that weights are properly normalized"""
        material = nbragg.materials["Fe_sg225_Iron-gamma.ncmat"]
        xs = nbragg.CrossSection.from_mtex(
            self.csv_file,
            material=material,
            short_name="gamma"
        )

        # Check that total weight equals 1
        total_weight = sum(mat['weight'] for mat in xs.materials.values())
        self.assertAlmostEqual(total_weight, 1.0, places=6)

        # Check individual weights
        # Original: 0.4 + 0.3 = 0.7, so powder should be 0.3
        self.assertAlmostEqual(xs.materials['gamma0']['weight'], 0.4, places=6)
        self.assertAlmostEqual(xs.materials['gamma1']['weight'], 0.3, places=6)
        self.assertAlmostEqual(xs.materials['gamma_powder']['weight'], 0.3, places=6)

    def test_from_mtex_orientation_parameters(self):
        """Test that orientation parameters are correctly set"""
        material = nbragg.materials["Fe_sg225_Iron-gamma.ncmat"]
        xs = nbragg.CrossSection.from_mtex(
            self.csv_file,
            material=material,
            short_name="gamma"
        )

        # Check first oriented phase
        phase0 = xs.materials['gamma0']

        # Check that dir1 and dir2 are lists of 3 elements
        self.assertIsInstance(phase0['dir1'], list)
        self.assertEqual(len(phase0['dir1']), 3)
        self.assertIsInstance(phase0['dir2'], list)
        self.assertEqual(len(phase0['dir2']), 3)

        # Check that directions are normalized (magnitude ≈ 1)
        import numpy as np
        dir1_mag = np.linalg.norm(phase0['dir1'])
        dir2_mag = np.linalg.norm(phase0['dir2'])
        self.assertAlmostEqual(dir1_mag, 1.0, places=6)
        self.assertAlmostEqual(dir2_mag, 1.0, places=6)

        # Check that mosaicity is set from fwhm
        self.assertAlmostEqual(phase0['mos'], 12.587233, places=6)

        # Check fixed parameters
        self.assertEqual(phase0['dirtol'], 1.0)
        self.assertEqual(phase0['theta'], 0.0)
        self.assertEqual(phase0['phi'], 0.0)

    def test_from_mtex_invalid_material_type(self):
        """Test that invalid material type raises TypeError"""
        with self.assertRaises(TypeError):
            nbragg.CrossSection.from_mtex(
                self.csv_file,
                material=123  # Invalid type
            )

    def test_from_mtex_nonexistent_file(self):
        """Test that nonexistent file raises FileNotFoundError"""
        material = nbragg.materials["Fe_sg225_Iron-gamma.ncmat"]
        with self.assertRaises(FileNotFoundError):
            nbragg.CrossSection.from_mtex(
                "nonexistent_file.csv",
                material=material
            )

    def test_from_mtex_material_properties_preserved(self):
        """Test that material properties like temp are preserved"""
        material = nbragg.materials["Fe_sg225_Iron-gamma.ncmat"].copy()
        material['temp'] = 350.0  # Custom temperature

        xs = nbragg.CrossSection.from_mtex(
            self.csv_file,
            material=material,
            short_name="gamma"
        )

        # Check that custom temperature is preserved
        self.assertEqual(xs.materials['gamma0']['temp'], 350.0)
        self.assertEqual(xs.materials['gamma1']['temp'], 350.0)
        self.assertEqual(xs.materials['gamma_powder']['temp'], 350.0)


if __name__ == '__main__':
    unittest.main()

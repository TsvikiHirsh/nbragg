"""Test suite for orientation parameter initialization in TransmissionModel"""
import unittest
import numpy as np
import nbragg


class TestOrientationParams(unittest.TestCase):
    """Test that orientation parameters are correctly initialized from materials"""

    def test_orientation_params_initialized_from_materials(self):
        """Test that θ, ϕ, η are initialized from material dictionary"""
        # Create a CrossSection with oriented materials
        materials = {
            'phase1': {
                'mat': 'Fe_sg229_Iron-alpha.ncmat',
                'temp': 300.0,
                'mos': 22.5,  # mosaicity in degrees
                'theta': 10.0,  # theta angle
                'phi': 15.0,  # phi angle
                'dir1': [1, 0, 0],
                'dir2': [0, 1, 0],
                'dirtol': 1.0,
                'weight': 0.5
            },
            'phase2': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'temp': 300.0,
                'mos': 18.3,
                'theta': 20.0,
                'phi': 25.0,
                'dir1': [0, 1, 0],
                'dir2': [0, 0, 1],
                'dirtol': 1.0,
                'weight': 0.5
            }
        }

        xs = nbragg.CrossSection(materials)

        # Create model with vary_orientation=False
        model = nbragg.TransmissionModel(
            xs,
            vary_orientation=False
        )

        # Check that orientation parameters exist and have correct initial values
        self.assertIn('θ_phase1', model.params)
        self.assertIn('ϕ_phase1', model.params)
        self.assertIn('η_phase1', model.params)
        self.assertIn('θ_phase2', model.params)
        self.assertIn('ϕ_phase2', model.params)
        self.assertIn('η_phase2', model.params)

        # Check that values are initialized from materials
        self.assertAlmostEqual(model.params['θ_phase1'].value, 10.0, places=6)
        self.assertAlmostEqual(model.params['ϕ_phase1'].value, 15.0, places=6)
        self.assertAlmostEqual(model.params['η_phase1'].value, 22.5, places=6)

        self.assertAlmostEqual(model.params['θ_phase2'].value, 20.0, places=6)
        self.assertAlmostEqual(model.params['ϕ_phase2'].value, 25.0, places=6)
        self.assertAlmostEqual(model.params['η_phase2'].value, 18.3, places=6)

        # Check that all are fixed (vary=False)
        self.assertFalse(model.params['θ_phase1'].vary)
        self.assertFalse(model.params['ϕ_phase1'].vary)
        self.assertFalse(model.params['η_phase1'].vary)
        self.assertFalse(model.params['θ_phase2'].vary)
        self.assertFalse(model.params['ϕ_phase2'].vary)
        self.assertFalse(model.params['η_phase2'].vary)

    def test_orientation_params_with_none_values(self):
        """Test that None values in materials default to 0"""
        materials = {
            'phase1': {
                'mat': 'Fe_sg229_Iron-alpha.ncmat',
                'temp': 300.0,
                'mos': None,  # Should default to 0
                'theta': None,  # Should default to 0
                'phi': None,  # Should default to 0
                'dir1': None,
                'dir2': None,
                'weight': 1.0
            }
        }

        xs = nbragg.CrossSection(materials)
        model = nbragg.TransmissionModel(xs, vary_orientation=False)

        # Check that None values defaulted to 0
        self.assertAlmostEqual(model.params['θ_phase1'].value, 0.0, places=6)
        self.assertAlmostEqual(model.params['ϕ_phase1'].value, 0.0, places=6)
        self.assertAlmostEqual(model.params['η_phase1'].value, 0.0, places=6)

    def test_orientation_params_with_missing_keys(self):
        """Test that missing keys in materials default to 0"""
        materials = {
            'phase1': {
                'mat': 'Fe_sg229_Iron-alpha.ncmat',
                'temp': 300.0,
                'weight': 1.0
                # No mos, theta, phi keys
            }
        }

        xs = nbragg.CrossSection(materials)
        model = nbragg.TransmissionModel(xs, vary_orientation=False)

        # Check that missing keys defaulted to 0
        self.assertAlmostEqual(model.params['θ_phase1'].value, 0.0, places=6)
        self.assertAlmostEqual(model.params['ϕ_phase1'].value, 0.0, places=6)
        self.assertAlmostEqual(model.params['η_phase1'].value, 0.0, places=6)

    def test_orientation_params_with_vary_true(self):
        """Test that vary=True is respected when vary_orientation=True"""
        materials = {
            'phase1': {
                'mat': 'Fe_sg229_Iron-alpha.ncmat',
                'temp': 300.0,
                'mos': 20.0,
                'theta': 5.0,
                'phi': 10.0,
                'dir1': [1, 0, 0],
                'dir2': [0, 1, 0],
                'weight': 1.0
            }
        }

        xs = nbragg.CrossSection(materials)
        model = nbragg.TransmissionModel(xs, vary_orientation=True)

        # Check that values are still initialized correctly
        self.assertAlmostEqual(model.params['θ_phase1'].value, 5.0, places=6)
        self.assertAlmostEqual(model.params['ϕ_phase1'].value, 10.0, places=6)
        self.assertAlmostEqual(model.params['η_phase1'].value, 20.0, places=6)

        # Check that vary=True
        self.assertTrue(model.params['θ_phase1'].vary)
        self.assertTrue(model.params['ϕ_phase1'].vary)
        self.assertTrue(model.params['η_phase1'].vary)

    def test_orientation_params_powder_phase(self):
        """Test orientation params for powder phases (should be 0 and fixed)"""
        materials = {
            'oriented': {
                'mat': 'Fe_sg229_Iron-alpha.ncmat',
                'temp': 300.0,
                'mos': 25.0,
                'theta': 15.0,
                'phi': 20.0,
                'dir1': [1, 0, 0],
                'dir2': [0, 1, 0],
                'weight': 0.7
            },
            'powder': {
                'mat': 'Fe_sg229_Iron-alpha.ncmat',
                'temp': 300.0,
                'mos': None,  # Powder has no orientation
                'theta': None,
                'phi': None,
                'dir1': None,
                'dir2': None,
                'weight': 0.3
            }
        }

        xs = nbragg.CrossSection(materials)
        model = nbragg.TransmissionModel(xs, vary_orientation=False)

        # Oriented phase should have non-zero values
        self.assertAlmostEqual(model.params['θ_oriented'].value, 15.0, places=6)
        self.assertAlmostEqual(model.params['ϕ_oriented'].value, 20.0, places=6)
        self.assertAlmostEqual(model.params['η_oriented'].value, 25.0, places=6)

        # Powder phase should have zero values
        self.assertAlmostEqual(model.params['θ_powder'].value, 0.0, places=6)
        self.assertAlmostEqual(model.params['ϕ_powder'].value, 0.0, places=6)
        self.assertAlmostEqual(model.params['η_powder'].value, 0.0, places=6)

    def test_orientation_params_from_mtex(self):
        """Test orientation params from MTEX-generated cross section"""
        import tempfile
        import pandas as pd
        import os

        # Create temporary MTEX CSV with orientation data
        temp_dir = tempfile.mkdtemp()
        csv_file = os.path.join(temp_dir, 'test_mtex.csv')

        test_data = {
            'alpha_mtex': [45.0, 90.0],
            'beta_mtex': [30.0, 45.0],
            'gamma_mtex': [60.0, 75.0],
            'volume_mtex': [0.3, 0.4],
            'xh': [1.0, 0.5],
            'xk': [0.0, 0.866],
            'xl': [0.0, 0.0],
            'yh': [0.0, -0.866],
            'yk': [1.0, 0.5],
            'yl': [0.0, 0.0],
            'fwhm': [15.5, 22.3]  # Mosaicity values
        }

        df = pd.DataFrame(test_data)
        df.to_csv(csv_file, index=False)

        try:
            # Create CrossSection from MTEX
            xs = nbragg.CrossSection.from_mtex(
                csv_file,
                material="Fe_sg225_Iron-gamma.ncmat",
                short_name="gamma",
                powder_phase=True
            )

            # Create model
            model = nbragg.TransmissionModel(xs, vary_orientation=False)

            # Check that mosaicity values were initialized correctly
            # from_mtex creates phases like gamma0, gamma1, gamma_powder
            self.assertIn('η_gamma0', model.params)
            self.assertIn('η_gamma1', model.params)

            # Check mosaicity values match MTEX fwhm
            self.assertAlmostEqual(model.params['η_gamma0'].value, 15.5, places=6)
            self.assertAlmostEqual(model.params['η_gamma1'].value, 22.3, places=6)

            # Theta and phi should be 0 (from_mtex sets them to 0)
            self.assertAlmostEqual(model.params['θ_gamma0'].value, 0.0, places=6)
            self.assertAlmostEqual(model.params['ϕ_gamma0'].value, 0.0, places=6)

            # Powder phase should have zero mosaicity
            self.assertAlmostEqual(model.params['η_gamma_powder'].value, 0.0, places=6)

        finally:
            # Clean up
            import shutil
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)


if __name__ == '__main__':
    unittest.main()

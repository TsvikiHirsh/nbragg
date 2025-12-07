"""
Tests for the new save/load API where result objects have .save() methods
and TransmissionModel.load() can handle both models and results.
"""

import unittest
import numpy as np
import pandas as pd
import tempfile
import os
import json
from pathlib import Path
from nbragg import CrossSection, TransmissionModel, materials


class TestNewSaveLoadAPI(unittest.TestCase):
    """Test the new save/load API."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.cross_section = CrossSection(iron=materials["Fe_sg229_Iron-alpha"])

        # Create mock data for fitting
        self.mock_data = pd.DataFrame({
            'wavelength': np.linspace(1, 5, 50),
            'trans': np.exp(-0.01 * np.sqrt(np.linspace(1, 5, 50))) + np.random.normal(0, 0.001, 50),
            'err': np.ones(50) * 0.01
        })

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    def test_result_has_save_method(self):
        """Test that fit results have a save() method."""
        model = TransmissionModel(self.cross_section, vary_background=True)
        result = model.fit(self.mock_data, wlmin=1.0, wlmax=5.0)

        # Check that result has a save method
        self.assertTrue(hasattr(result, 'save'))
        self.assertTrue(callable(result.save))

    def test_result_save_creates_files(self):
        """Test that result.save() creates both result and model files."""
        model = TransmissionModel(self.cross_section, vary_background=True)
        result = model.fit(self.mock_data, wlmin=1.0, wlmax=5.0)

        result_file = os.path.join(self.temp_dir, 'test_result.json')
        model_file = os.path.join(self.temp_dir, 'test_result_model.json')

        result.save(result_file)

        # Check both files were created
        self.assertTrue(os.path.exists(result_file))
        self.assertTrue(os.path.exists(model_file))

    def test_model_load_from_model_file(self):
        """Test loading a model from a model file."""
        model = TransmissionModel(self.cross_section, vary_background=True, tof_length=12.5)
        model_file = os.path.join(self.temp_dir, 'test_model.json')

        model.save(model_file)
        loaded_model = TransmissionModel.load(model_file)

        # Check model attributes
        self.assertEqual(loaded_model.tof_length, 12.5)
        self.assertIsNotNone(loaded_model.background)

    def test_model_load_from_result_file(self):
        """Test loading a model from a result file."""
        model = TransmissionModel(self.cross_section, vary_background=True)
        result = model.fit(self.mock_data, wlmin=1.0, wlmax=5.0)

        result_file = os.path.join(self.temp_dir, 'test_result.json')
        result.save(result_file)

        loaded_model = TransmissionModel.load(result_file)

        # Check that model has a result attribute
        self.assertTrue(hasattr(loaded_model, 'result'))
        self.assertIsNotNone(loaded_model.result)

        # Check result attributes
        self.assertAlmostEqual(loaded_model.result.redchi, result.redchi, places=10)
        self.assertAlmostEqual(loaded_model.result.chisqr, result.chisqr, places=10)

    def test_loaded_result_has_all_methods(self):
        """Test that loaded result has plot methods and other attributes."""
        model = TransmissionModel(self.cross_section, vary_background=True)
        result = model.fit(self.mock_data, wlmin=1.0, wlmax=5.0)

        result_file = os.path.join(self.temp_dir, 'test_result.json')
        result.save(result_file)

        loaded_model = TransmissionModel.load(result_file)

        # Check that loaded result has the expected methods
        self.assertTrue(hasattr(loaded_model.result, 'plot'))
        self.assertTrue(hasattr(loaded_model.result, 'plot_total_xs'))
        self.assertTrue(hasattr(loaded_model.result, 'show_available_params'))
        self.assertTrue(hasattr(loaded_model.result, 'save'))
        self.assertTrue(callable(loaded_model.result.save))

    def test_loaded_result_parameters_match(self):
        """Test that loaded result has the correct fitted parameters."""
        model = TransmissionModel(self.cross_section, vary_background=True)
        result = model.fit(self.mock_data, wlmin=1.0, wlmax=5.0)

        result_file = os.path.join(self.temp_dir, 'test_result.json')
        result.save(result_file)

        loaded_model = TransmissionModel.load(result_file)

        # Check parameters match
        for key in result.params.keys():
            self.assertAlmostEqual(
                result.params[key].value,
                loaded_model.result.params[key].value,
                places=10
            )

    def test_init_from_model_file(self):
        """Test initializing TransmissionModel from a model file."""
        model = TransmissionModel(self.cross_section, vary_background=True, tof_length=15.0)
        model_file = os.path.join(self.temp_dir, 'test_model.json')

        model.save(model_file)
        loaded_model = TransmissionModel(model_file)

        self.assertEqual(loaded_model.tof_length, 15.0)
        self.assertIsNotNone(loaded_model.background)

    def test_init_from_result_file(self):
        """Test initializing TransmissionModel from a result file."""
        model = TransmissionModel(self.cross_section, vary_background=True)
        result = model.fit(self.mock_data, wlmin=1.0, wlmax=5.0)

        result_file = os.path.join(self.temp_dir, 'test_result.json')
        result.save(result_file)

        loaded_model = TransmissionModel(result_file)

        # Check that model has a result attribute
        self.assertTrue(hasattr(loaded_model, 'result'))
        self.assertIsNotNone(loaded_model.result)
        self.assertAlmostEqual(loaded_model.result.redchi, result.redchi, places=10)

    def test_stages_preserved_in_model(self):
        """Test that stages are preserved when saving/loading a model."""
        model = TransmissionModel(self.cross_section, vary_background=True, vary_response=True)
        model.stages = {'stage1': ['norm', 'thickness'], 'stage2': 'background'}

        model_file = os.path.join(self.temp_dir, 'test_model.json')
        model.save(model_file)

        loaded_model = TransmissionModel.load(model_file)

        self.assertEqual(loaded_model.stages, model.stages)

    def test_loaded_result_can_be_saved_again(self):
        """Test that a loaded result can be saved again."""
        model = TransmissionModel(self.cross_section, vary_background=True)
        result = model.fit(self.mock_data, wlmin=1.0, wlmax=5.0)

        result_file1 = os.path.join(self.temp_dir, 'result1.json')
        result.save(result_file1)

        loaded_model = TransmissionModel.load(result_file1)

        result_file2 = os.path.join(self.temp_dir, 'result2.json')
        loaded_model.result.save(result_file2)

        # Both result files should exist
        self.assertTrue(os.path.exists(result_file1))
        self.assertTrue(os.path.exists(result_file2))

        # Load the second result and check it matches
        loaded_model2 = TransmissionModel.load(result_file2)
        self.assertAlmostEqual(
            loaded_model.result.redchi,
            loaded_model2.result.redchi,
            places=10
        )

    def test_multiphase_save_load(self):
        """Test save/load with multiphase materials."""
        xs = CrossSection({
            'alpha': {
                'mat': 'Fe_sg229_Iron-alpha.ncmat',
                'weight': 0.4
            },
            'gamma': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                'weight': 0.6
            }
        })

        model = TransmissionModel(xs, vary_background=True)
        result = model.fit(self.mock_data, wlmin=1.0, wlmax=5.0)

        result_file = os.path.join(self.temp_dir, 'multiphase_result.json')
        result.save(result_file)

        loaded_model = TransmissionModel.load(result_file)

        # Check materials are preserved
        self.assertEqual(len(loaded_model._materials), 2)
        self.assertIn('alpha', loaded_model._materials)
        self.assertIn('gamma', loaded_model._materials)

        # Check result is preserved
        self.assertIsNotNone(loaded_model.result)
        self.assertAlmostEqual(loaded_model.result.redchi, result.redchi, places=10)

    def test_model_from_result_can_fit_again(self):
        """Test that a model loaded from a result can be used to fit again."""
        model = TransmissionModel(self.cross_section, vary_background=True)
        result1 = model.fit(self.mock_data, wlmin=1.0, wlmax=5.0)

        result_file = os.path.join(self.temp_dir, 'result.json')
        result1.save(result_file)

        loaded_model = TransmissionModel.load(result_file)

        # Fit again with the loaded model
        result2 = loaded_model.fit(self.mock_data, wlmin=1.0, wlmax=5.0)

        # Should be able to fit successfully
        self.assertTrue(result2.success)
        self.assertIsNotNone(result2.redchi)

    def test_result_file_format(self):
        """Test that result files have the correct format."""
        model = TransmissionModel(self.cross_section, vary_background=True)
        result = model.fit(self.mock_data, wlmin=1.0, wlmax=5.0)

        result_file = os.path.join(self.temp_dir, 'result.json')
        result.save(result_file)

        with open(result_file, 'r') as f:
            data = json.load(f)

        self.assertEqual(data['version'], '1.0')
        self.assertEqual(data['class'], 'ModelResult')
        self.assertIn('params', data)
        self.assertIn('chisqr', data)
        self.assertIn('redchi', data)
        self.assertIn('success', data)

    def test_model_file_format(self):
        """Test that model files have the correct format."""
        model = TransmissionModel(self.cross_section, vary_background=True, tof_length=10.0)

        model_file = os.path.join(self.temp_dir, 'model.json')
        model.save(model_file)

        with open(model_file, 'r') as f:
            data = json.load(f)

        self.assertEqual(data['version'], '1.0')
        self.assertEqual(data['class'], 'TransmissionModel')
        self.assertIn('materials', data)
        self.assertIn('params', data)
        self.assertIn('stages', data)
        self.assertIn('tof_length', data)
        self.assertEqual(data['tof_length'], 10.0)

    def test_loaded_result_has_required_attributes(self):
        """Test that loaded result has all attributes required by lmfit."""
        model = TransmissionModel(self.cross_section, vary_background=True)
        result = model.fit(self.mock_data, wlmin=1.0, wlmax=5.0)

        result_file = os.path.join(self.temp_dir, 'result.json')
        result.save(result_file)

        loaded_model = TransmissionModel.load(result_file)

        # Check for attributes that lmfit's _repr_html_ needs
        self.assertTrue(hasattr(loaded_model.result, 'method'))
        self.assertTrue(hasattr(loaded_model.result, 'var_names'))
        self.assertTrue(hasattr(loaded_model.result, 'init_vals'))
        self.assertTrue(hasattr(loaded_model.result, 'aborted'))
        self.assertTrue(hasattr(loaded_model.result, 'errorbars'))

        # Try to call _repr_html_ (would fail with AttributeError before fix)
        html = loaded_model.result._repr_html_()
        self.assertIsInstance(html, str)
        self.assertGreater(len(html), 0)

    def test_loaded_model_has_fitted_values(self):
        """Test that loaded model uses fitted parameter values, not initial values."""
        model = TransmissionModel(self.cross_section, vary_background=True)
        result = model.fit(self.mock_data, wlmin=1.0, wlmax=5.0)

        # Store fitted values
        fitted_thickness = result.params['thickness'].value
        fitted_norm = result.params['norm'].value

        result_file = os.path.join(self.temp_dir, 'result.json')
        result.save(result_file)

        loaded_model = TransmissionModel.load(result_file)

        # Loaded model should have fitted values in both model.params and model.result.params
        self.assertAlmostEqual(loaded_model.params['thickness'].value, fitted_thickness, places=10)
        self.assertAlmostEqual(loaded_model.params['norm'].value, fitted_norm, places=10)
        self.assertAlmostEqual(loaded_model.result.params['thickness'].value, fitted_thickness, places=10)
        self.assertAlmostEqual(loaded_model.result.params['norm'].value, fitted_norm, places=10)


if __name__ == '__main__':
    unittest.main()

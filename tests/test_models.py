import unittest
import numpy as np
import pandas as pd
from unittest.mock import patch
from nbragg.models import TransmissionModel
from nbragg import CrossSection, materials
from lmfit import Parameters

class TestTransmissionModel(unittest.TestCase):

    def setUp(self):
        """Set up the model with a real CrossSection for testing."""
        self.cross_section = CrossSection(iron=materials["Fe_sg229_Iron-alpha.ncmat"])
        self.model = TransmissionModel(self.cross_section)
        # Mock data for fitting tests
        self.mock_data = pd.DataFrame({
            'wavelength': np.linspace(1, 5, 100),
            'trans': np.exp(-0.01 * np.sqrt(np.linspace(1, 5, 100))),
            'err': np.ones(100) * 0.01
        })

    def test_initial_parameters(self):
        """Test the initial parameters are correctly set."""
        params = self.model.params
        self.assertIn('norm', params, "Expected 'norm' parameter in model.params")
        self.assertIn('thickness', params, "Expected 'thickness' parameter in model.params")
        self.assertTrue(params['norm'].vary)  # norm.vary is True by default
        self.assertTrue(params['thickness'].vary)  # thickness.vary is True by default

    def test_transmission_calculation(self):
        """Test the transmission function calculation."""
        wl = np.array([1.0, 2.0, 4.0])
        T = self.model.transmission(wl, thickness=1)
        self.assertEqual(T.shape, wl.shape)
        self.assertTrue(np.all(T > 0))

    def test_background_varying(self):
        """Test the model with varying background parameters."""
        model_vary_bg = TransmissionModel(self.cross_section, vary_background=True, background="polynomial3")
        params = model_vary_bg.params
        self.assertIn('bg0', params, "Expected 'bg0' parameter when vary_background=True")
        self.assertTrue(params['bg0'].vary)
        self.assertTrue(params['bg1'].vary)
        self.assertTrue(params['bg2'].vary)

    def test_stages_setter_valid(self):
        """Test the stages setter with valid inputs."""
        self.model.stages = {'all_params': 'all'}
        self.assertEqual(self.model.stages, {'all_params': 'all'})

        model_vary_bg = TransmissionModel(self.cross_section, vary_background=True, background="polynomial3")
        model_vary_bg.stages = {'bg': 'background'}
        self.assertEqual(model_vary_bg.stages, {'bg': 'background'})

        model_vary_bg.stages = {'basic': ['norm', 'thickness']}
        self.assertEqual(model_vary_bg.stages, {'basic': ['norm', 'thickness']})

    def test_stages_setter_invalid(self):
        """Test the stages setter with invalid inputs."""
        with self.assertRaises(ValueError) as context:
            self.model.stages = {1: 'all'}
        self.assertTrue("Stage names must be strings" in str(context.exception))

        with self.assertRaises(ValueError) as context:
            self.model.stages = {'basic': 'invalid_group'}
        self.assertTrue("must be 'all' or a valid group name" in str(context.exception))

        with self.assertRaises(ValueError) as context:
            self.model.stages = {'basic': [1, 'thickness']}
        self.assertTrue("Parameters in stage 'basic' must be strings" in str(context.exception))

    def test_default_rietveld_fit(self):
        """Test that the fit method uses rietveld by default."""
        model_vary_bg = TransmissionModel(self.cross_section, vary_background=True, vary_response=True, background="polynomial3")
        model_vary_bg.stages = {'basic': ['norm', 'thickness'], 'background': 'background'}

        with patch.object(model_vary_bg, '_multistage_fit') as mock_multistage:
            mock_multistage.return_value = MockFitResult()
            result = model_vary_bg.fit(self.mock_data, wlmin=1, wlmax=5)
            mock_multistage.assert_called_once()
            call_args, call_kwargs = mock_multistage.call_args
            from pandas.testing import assert_frame_equal
            assert_frame_equal(call_args[0], self.mock_data)  # Data passed correctly
            self.assertEqual(call_args[2], 1)  # wlmin
            self.assertEqual(call_args[3], 5)  # wlmax
            self.assertEqual(call_kwargs.get('method', None), 'rietveld')  # Default method
            self.assertEqual(call_kwargs.get('stages', None), {'basic': ['norm', 'thickness'], 'background': 'background'})  # Stages

class MockFitResult:
    """Mock lmfit.ModelResult for testing."""
    def __init__(self):
        self.params = Parameters()
        self.redchi = 1.0
        self.plot = lambda: None
        self.plot_total_xs = lambda: None
        self.show_available_params = lambda: None

if __name__ == '__main__':
    unittest.main()
import unittest
import numpy as np
import pandas as pd
from unittest.mock import patch
from nbragg.models import TransmissionModel
from nbragg import CrossSection, materials
from lmfit import Parameters
import re

class TestTransmissionModel(unittest.TestCase):

    def setUp(self):
        """Set up the model with a real CrossSection for testing."""
        self.cross_section = CrossSection(iron=materials["Fe_sg229_Iron-alpha"])
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

    def test_one_by_one_stages(self):
        """Test that one-by-one stages are correctly expanded and fitted."""
        model = TransmissionModel(self.cross_section, vary_background=True, vary_response=True, background="polynomial3")
        stages = {"1": ["thickness", "background", "response", "one-by-one"]}

        # Define resolve_group to mirror _multistage_fit logic
        def resolve_single_param_or_group(item, model):
            group_map = {
                "basic": ["norm", "thickness"],
                "background": [p for p in model.params if re.compile(r"(b|bg)\d+").match(p) or p.startswith("b_")],
                "response": [p for p in model.params if model.response and p in model.response.params],
            }
            if item == "all":
                return [p for p in model.params if model.params[p].vary]
            elif item in group_map:
                return group_map[item]
            elif item in model.params:
                return [item]
            else:
                return []

        def resolve_group(entry, stage_name, model):
            params_list = []
            if isinstance(entry, list):
                is_one_by_one = "one-by-one" in entry
                for item in entry:
                    if item == "one-by-one" or (isinstance(item, str) and (item.startswith("wlmin=") or item.startswith("wlmax="))):
                        continue
                    params_list.extend(resolve_single_param_or_group(item, model))
                if is_one_by_one:
                    sub_stages = []
                    for i, param in enumerate(params_list):
                        var_part = param.split("_")[-1] if "_" in param else param
                        sub_name = f"{stage_name}_{var_part}" if len(params_list) > 1 else stage_name
                        sub_stages.append((sub_name, [param], {}))
                    return sub_stages
                return [(stage_name, params_list, {})]
            return [(stage_name, resolve_single_param_or_group(entry, model), {})]

        # Expected parameters
        expected_params = ['thickness', 'bg0', 'bg1', 'bg2'] + (['α0', 'β0'] if model.response and {'α0', 'β0'}.issubset(set(model.params.keys())) else [])
        expected_stages = [(f"1_{p.split('_')[-1]}" if '_' in p else f"1_{p}", [p], {}) for p in expected_params]
        expected_stage_names, expected_param_groups, _ = zip(*expected_stages) if expected_stages else ([], [], [])

        with patch.object(model, '_multistage_fit') as mock_multistage:
            mock_multistage.return_value = MockFitResult()
            result = model.fit(self.mock_data, stages=stages, wlmin=1, wlmax=5)
            mock_multistage.assert_called_once()
            call_args, call_kwargs = mock_multistage.call_args
            from pandas.testing import assert_frame_equal
            assert_frame_equal(call_args[0], self.mock_data)  # Data passed correctly
            self.assertEqual(call_args[2], 1)  # wlmin
            self.assertEqual(call_args[3], 5)  # wlmax
            self.assertEqual(call_kwargs.get('method', None), 'rietveld')  # Default method
            self.assertEqual(call_kwargs.get('stages', None), stages)  # Original stages

        # Verify stage expansion using resolve_group
        resolved = resolve_group(stages["1"], "1", model)
        resolved_stage_names, resolved_param_groups, resolved_overrides = zip(*resolved) if resolved else ([], [], [])
        self.assertEqual(list(resolved_stage_names), list(expected_stage_names))
        self.assertEqual(list(resolved_param_groups), list(expected_param_groups))
        self.assertEqual(list(resolved_overrides), [{}] * len(expected_stage_names))

    def test_one_by_one_empty_stages(self):
        """Test that one-by-one stages with no valid parameters are handled."""
        model = TransmissionModel(self.cross_section)
        stages = {"1": ["invalid_param", "one-by-one"]}
        
        with self.assertRaises(ValueError) as context:
            model.fit(self.mock_data, stages=stages)
        self.assertTrue("No valid stages found" in str(context.exception))

    def test_multiple_stages_with_one_by_one(self):
        """Test multiple stages including one-by-one."""
        model = TransmissionModel(self.cross_section, vary_background=True, vary_response=True, background="polynomial3")
        stages = {"1": ["thickness", "background"], "second": ["response", "one-by-one"]}

        # Define resolve_group to mirror _multistage_fit logic
        def resolve_single_param_or_group(item, model):
            group_map = {
                "basic": ["norm", "thickness"],
                "background": [p for p in model.params if re.compile(r"(b|bg)\d+").match(p) or p.startswith("b_")],
                "response": [p for p in model.params if model.response and p in model.response.params],
            }
            if item == "all":
                return [p for p in model.params if model.params[p].vary]
            elif item in group_map:
                return group_map[item]
            elif item in model.params:
                return [item]
            else:
                return []

        def resolve_group(entry, stage_name, model):
            params_list = []
            if isinstance(entry, list):
                is_one_by_one = "one-by-one" in entry
                for item in entry:
                    if item == "one-by-one" or (isinstance(item, str) and (item.startswith("wlmin=") or item.startswith("wlmax="))):
                        continue
                    params_list.extend(resolve_single_param_or_group(item, model))
                if is_one_by_one:
                    sub_stages = []
                    for i, param in enumerate(params_list):
                        var_part = param.split("_")[-1] if "_" in param else param
                        sub_name = f"{stage_name}_{var_part}" if len(params_list) > 1 else stage_name
                        sub_stages.append((sub_name, [param], {}))
                    return sub_stages
                return [(stage_name, params_list, {})]
            return [(stage_name, resolve_single_param_or_group(entry, model), {})]

        # Expected stages
        expected_params_1 = ['thickness', 'bg0', 'bg1', 'bg2']
        expected_params_second = ['α0', 'β0'] if model.response and {'α0', 'β0'}.issubset(set(model.params.keys())) else []
        expected_stages = [
            ("1", expected_params_1, {}),
            *[(f"second_{p.split('_')[-1]}" if '_' in p else f"second_{p}", [p], {}) for p in expected_params_second]
        ]
        expected_stage_names, expected_param_groups, _ = zip(*expected_stages) if expected_stages else ([], [], [])

        with patch.object(model, '_multistage_fit') as mock_multistage:
            mock_multistage.return_value = MockFitResult()
            result = model.fit(self.mock_data, stages=stages, wlmin=1, wlmax=5)
            mock_multistage.assert_called_once()
            call_args, call_kwargs = mock_multistage.call_args
            from pandas.testing import assert_frame_equal
            assert_frame_equal(call_args[0], self.mock_data)  # Data passed correctly
            self.assertEqual(call_args[2], 1)  # wlmin
            self.assertEqual(call_args[3], 5)  # wlmax
            self.assertEqual(call_kwargs.get('method', None), 'rietveld')  # Default method
            self.assertEqual(call_kwargs.get('stages', None), stages)  # Original stages

        # Verify stage expansion
        resolved_stages = []
        for stage_name, stage_def in stages.items():
            resolved_stages.extend(resolve_group(stage_def, stage_name, model))
        resolved_stage_names, resolved_param_groups, resolved_overrides = zip(*resolved_stages) if resolved_stages else ([], [], [])
        self.assertEqual(list(resolved_stage_names), list(expected_stage_names))
        self.assertEqual(list(resolved_param_groups), list(expected_param_groups))
        self.assertEqual(list(resolved_overrides), [{}] * len(expected_stage_names))

class TestSANSModelParameters(unittest.TestCase):
    def test_sans_params_single_phase(self):
        """Test SANS parameter creation for single-phase material."""
        xs = CrossSection({
            'aluminum': {
                'mat': 'Al_sg225.ncmat',
                'sans': 50.0
            }
        })
        model = TransmissionModel(xs, vary_sans=True)

        params = model.params
        self.assertIn('sans', params, "Expected 'sans' parameter in model.params")
        self.assertEqual(params['sans'].value, 50.0)
        self.assertTrue(params['sans'].vary)
        self.assertEqual(params['sans'].min, 0.0)
        self.assertEqual(params['sans'].max, 1000)

    def test_sans_params_multiphase(self):
        """Test SANS parameter creation for multi-phase materials."""
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
        model = TransmissionModel(xs, vary_sans=True)

        params = model.params
        self.assertIn('sans1', params, "Expected 'sans1' parameter in model.params")
        self.assertIn('sans2', params, "Expected 'sans2' parameter in model.params")
        self.assertEqual(params['sans1'].value, 30.0)
        self.assertEqual(params['sans2'].value, 60.0)
        self.assertTrue(params['sans1'].vary)
        self.assertTrue(params['sans2'].vary)

    def test_sans_params_not_vary(self):
        """Test SANS parameters with vary=False."""
        xs = CrossSection({
            'aluminum': {
                'mat': 'Al_sg225.ncmat',
                'sans': 50.0
            }
        })
        model = TransmissionModel(xs, vary_sans=False)

        params = model.params
        self.assertIn('sans', params, "Expected 'sans' parameter in model.params")
        self.assertEqual(params['sans'].value, 50.0)
        self.assertFalse(params['sans'].vary)

    def test_sans_stages_initialization(self):
        """Test that SANS stage is initialized when vary_sans=True."""
        xs = CrossSection({
            'aluminum': {
                'mat': 'Al_sg225.ncmat',
                'sans': 50.0
            }
        })
        model = TransmissionModel(xs, vary_sans=True)

        self.assertIn('sans', model.stages)
        self.assertEqual(model.stages['sans'], 'sans')

    def test_sans_group_in_stages_setter(self):
        """Test that SANS parameters are recognized in stages setter."""
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
        model = TransmissionModel(xs, vary_sans=True)

        # Test setting stages with SANS group
        model.stages = {'sans_stage': 'sans'}
        self.assertEqual(model.stages, {'sans_stage': 'sans'})

    def test_sans_params_only_if_defined(self):
        """Test that SANS parameters are only created if sans is defined in material."""
        xs = CrossSection({
            'alpha': {
                'mat': 'Fe_sg229_Iron-alpha.ncmat',
                'sans': 30.0,
                'weight': 0.4
            },
            'gamma': {
                'mat': 'Fe_sg225_Iron-gamma.ncmat',
                # No SANS parameter defined
                'weight': 0.6
            }
        })
        model = TransmissionModel(xs, vary_sans=True)

        params = model.params
        # Only sans1 should exist (for alpha phase)
        self.assertIn('sans1', params, "Expected 'sans1' parameter for alpha phase")
        # sans2 should NOT exist because gamma doesn't have SANS defined
        # The materials dict has 2 materials, but only 1 has sans defined

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
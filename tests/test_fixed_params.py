"""
Test that parameters can be fixed (vary=False) during fitting.
"""
import pytest
import numpy as np
import pandas as pd
import nbragg


def test_fix_norm_parameter():
    """Test that setting norm.vary=False prevents it from being fitted."""
    # Create simple mock data
    wavelength = np.linspace(1.0, 5.0, 100)
    trans = 0.9 * np.exp(-0.1 * wavelength)  # norm should be ~0.9
    err = 0.01 * np.ones_like(trans)

    data = pd.DataFrame({
        'wavelength': wavelength,
        'trans': trans,
        'err': err
    })

    # Create a simple model
    xs = nbragg.CrossSection(iron="Fe_sg229_Iron-alpha.ncmat")
    model = nbragg.TransmissionModel(xs, vary_basic=True)

    # Set initial value for norm and fix it
    initial_norm = 0.85
    model.params["norm"].set(value=initial_norm, vary=False)

    # Verify it's set correctly before fit
    assert model.params["norm"].value == initial_norm
    assert model.params["norm"].vary == False

    # Fit the model
    result = model.fit(data, wlmin=1.0, wlmax=5.0)

    # Check that norm did NOT change
    assert result.params["norm"].value == initial_norm, \
        f"norm changed from {initial_norm} to {result.params['norm'].value}, but vary=False!"

    # Check that thickness DID vary (it should be in the fit)
    assert result.params["thickness"].vary == True
    assert result.params["thickness"].value != model.params["thickness"].init_value


def test_fix_thickness_parameter():
    """Test that setting thickness.vary=False prevents it from being fitted."""
    # Create simple mock data
    wavelength = np.linspace(1.0, 5.0, 100)
    trans = np.exp(-0.1 * wavelength)
    err = 0.01 * np.ones_like(trans)

    data = pd.DataFrame({
        'wavelength': wavelength,
        'trans': trans,
        'err': err
    })

    # Create a simple model
    xs = nbragg.CrossSection(iron="Fe_sg229_Iron-alpha.ncmat")
    model = nbragg.TransmissionModel(xs, vary_basic=True)

    # Set initial value for thickness and fix it
    initial_thickness = 0.05
    model.params["thickness"].set(value=initial_thickness, vary=False)

    # Verify it's set correctly before fit
    assert model.params["thickness"].value == initial_thickness
    assert model.params["thickness"].vary == False

    # Fit the model
    result = model.fit(data, wlmin=1.0, wlmax=5.0)

    # Check that thickness did NOT change
    assert result.params["thickness"].value == initial_thickness, \
        f"thickness changed from {initial_thickness} to {result.params['thickness'].value}, but vary=False!"

    # Check that norm DID vary (it should be in the fit)
    assert result.params["norm"].vary == True


def test_fix_multiple_parameters():
    """Test that multiple parameters can be fixed simultaneously."""
    # Create simple mock data
    wavelength = np.linspace(1.0, 5.0, 100)
    trans = 0.9 * np.exp(-0.1 * wavelength)
    err = 0.01 * np.ones_like(trans)

    data = pd.DataFrame({
        'wavelength': wavelength,
        'trans': trans,
        'err': err
    })

    # Create a model with background
    xs = nbragg.CrossSection(iron="Fe_sg229_Iron-alpha.ncmat")
    model = nbragg.TransmissionModel(xs, vary_basic=True, vary_background=True)

    # Fix norm and one background parameter
    initial_norm = 0.85
    initial_bg0 = 0.05
    model.params["norm"].set(value=initial_norm, vary=False)
    model.params["bg0"].set(value=initial_bg0, vary=False)

    # Fit the model
    result = model.fit(data, wlmin=1.0, wlmax=5.0)

    # Check that fixed parameters did NOT change
    assert result.params["norm"].value == initial_norm
    assert result.params["bg0"].value == initial_bg0

    # Check that other parameters DID vary
    assert result.params["thickness"].vary == True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

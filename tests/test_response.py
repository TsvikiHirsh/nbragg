"""
Tests for the Jorgensen peak-profile implementations in nbragg.response.

The reference is a verbatim copy of `Jorgensen_simple` from ORNL's
braggedgemodeling (bem/peak_profile.py), which any working Jorgensen
implementation must match up to normalization.
"""
import numpy as np
import pytest
from scipy.special import erfc

from nbragg.response import Response


def jorgensen_simple_ref(x, sigma, alpha, beta):
    """Verbatim ORNL bem.peak_profile.Jorgensen_simple (un-normalized)."""
    scale = alpha * beta / 2 / (alpha + beta)
    sigma2 = sigma * sigma
    sqrt2 = np.sqrt(2)
    u = alpha / 2. * (alpha * sigma2 + 2 * x)
    v = beta / 2. * (beta * sigma2 - 2 * x)
    y = (alpha * sigma2 + x) / (sqrt2 * sigma)
    z = (beta * sigma2 - x) / (sqrt2 * sigma)
    with np.errstate(over="ignore", invalid="ignore"):
        term1 = np.exp(u) * erfc(y)
        term1[erfc(y) == 0] = 0
        term2 = np.exp(v) * erfc(z)
        term2[erfc(z) == 0] = 0
    return scale * (term1 + term2)


def _normalize(p):
    s = p.sum()
    return p / s if s > 0 else p


class TestFullJorgensenAgainstReference:
    """full_jorgensen_response must match ORNL Jorgensen_simple verbatim."""

    @pytest.mark.parametrize(
        "wl, α0, α1, β0, β1, σ0, σ1, σ2",
        [
            # Friend's converged values (TOF imaging, GSAS-style)
            (4.0, 0.0, 494.32, 48.6, 0.0, 0.0, 1.0e-5, 0.0),
            # σ1 only (response-only)
            (4.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0e-4, 0.0),
            # Mixed α0/α1
            (3.0, 1.0, 100.0, 50.0, 0.0, 0.0, 5.0e-5, 0.0),
            # σ0 + σ1 + σ2 all on
            (5.0, 0.0, 200.0, 30.0, 5.0, 1e-6, 2e-5, 1e-6),
        ],
    )
    def test_matches_ornl_reference(self, wl, α0, α1, β0, β1, σ0, σ1, σ2):
        r = Response(kind="full_jorgensen")
        ours = r.function(wl=wl, α0=α0, α1=α1, β0=β0, β1=β1, σ0=σ0, σ1=σ1, σ2=σ2)

        # Recompute the reference on the same grid with the same d-dependence
        d = wl / 2.0
        alpha = α0 + α1 / d
        beta = β0 + β1 / (d ** 4)
        sigma = np.sqrt(σ0 ** 2 + (σ1 * d) ** 2 + (σ2 * d * d) ** 2)
        ref = jorgensen_simple_ref(r.tgrid, sigma, alpha, beta)
        ref = np.nan_to_num(ref, nan=0.0, posinf=0.0, neginf=0.0)
        ref = _normalize(ref)

        np.testing.assert_allclose(ours, ref, rtol=1e-10, atol=1e-12)


class TestFullJorgensenBasicProperties:
    def setup_method(self):
        self.r = Response(kind="full_jorgensen")

    def test_normalized(self):
        p = self.r.function(wl=4.0, **self.r.params.valuesdict())
        assert p.sum() == pytest.approx(1.0, rel=1e-10)

    def test_non_negative(self):
        p = self.r.function(wl=4.0, **self.r.params.valuesdict())
        assert (p >= 0).all()

    def test_unimodal(self):
        p = self.r.function(wl=4.0, **self.r.params.valuesdict())
        # Exactly one sign change from + to - in the first derivative
        sign_changes = (np.diff(np.sign(np.diff(p))) < 0).sum()
        assert sign_changes == 1

    def test_sigma1_widens_profile(self):
        """σ1 should be a real knob: increasing it makes the profile wider."""
        base = self.r.params.valuesdict()
        widths = []
        for s1 in [5e-6, 1e-5, 5e-5, 1e-4]:
            kw = {**base, "σ1": s1}
            p = self.r.function(wl=4.0, **kw)
            widths.append((p > 0.5 * p.max()).sum())
        # Monotone non-decreasing — σ1 must actually do something
        for a, b in zip(widths, widths[1:]):
            assert b >= a, f"FWHM samples went down: {widths}"
        # At least one of the steps must strictly widen the profile
        assert widths[-1] > widths[0], f"σ1 had no measurable effect: {widths}"

    def test_alpha_plus_beta_must_be_positive(self):
        with pytest.raises(ValueError):
            self.r.function(wl=4.0, α0=-100.0, α1=0.0, β0=10.0, β1=0.0,
                            σ0=0.0, σ1=1e-5, σ2=0.0)

    def test_staged_preset_shape(self):
        preset = self.r.staged_preset()
        assert list(preset.keys()) == ["response_sig1", "response_beta", "response_alpha"]
        assert preset["response_sig1"] == ["σ1"]
        assert "σ1" in preset["response_beta"] and "β0" in preset["response_beta"]
        assert "α1" in preset["response_alpha"]


class TestSimpleJorgensenStillWorks:
    """Smoke test: the simple `jorgensen` kind isn't broken by these changes."""

    def test_simple_jorgensen_normalized(self):
        r = Response(kind="jorgensen")
        p = r.function(**r.params.valuesdict())
        assert p.sum() == pytest.approx(1.0, rel=1e-10)
        assert (p >= 0).all()

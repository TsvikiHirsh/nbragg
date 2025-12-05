# plot_total_xs Enhancements

## Summary

Enhanced the `plot_total_xs` method to display actual data and individual phase contributions, making it much more useful for analyzing multi-phase materials.

## Issues Fixed

### 1. Data Not Visible in Plot

**Problem:** The `plot_total_xs` method only showed the model curve, not the actual experimental data points.

**Solution:**
- Extract transmission data from `fit_result.data` and `fit_result.weights`
- Convert transmission back to cross section using the inverse relationship
- Plot data with error bars in light purple (`#B0A1BA`)
- Error bars properly calculated using error propagation

### 2. No Way to Show Individual Phase Contributions

**Problem:** For multi-phase materials, there was no way to visualize how each phase contributes to the total cross section.

**Solution:** Added `split_phases=True` option that:
- Plots each phase with a different color from the `turbo` colormap
- Shows weighted contributions: `xs_phase * weight_phase`
- Adds text labels on the right showing phase percentages
- Sorts phases by weight (largest first)
- Filters out negligible contributions (< 1% of total)

### 3. No Residuals Visualization

**Problem:** Users had no quick way to assess fit quality in cross section space.

**Solution:** Added `plot_residuals=True` option that:
- Creates a 2-panel plot (cross section + residuals)
- Bottom panel shows `data_xs - model_xs`
- Includes a zero reference line
- Properly shares x-axis between panels

## New Features

### 1. Data Plotting

When a fit has been performed, data points are automatically shown:

```python
result.plot_total_xs()  # Now shows data + model
```

Features:
- Error bars displayed
- Legend shows chi-squared value
- Data appears behind model line for clarity

### 2. Split Phases Visualization

View individual phase contributions:

```python
result.plot_total_xs(split_phases=True)
```

Features:
- Each phase plotted with unique color
- Phase weights shown as percentages
- Labels positioned at right edge of plot
- Automatic sorting by weight

### 3. Residuals Panel

Assess fit quality:

```python
ax_main, ax_res = result.plot_total_xs(plot_residuals=True)
```

Features:
- Two-panel layout (4:1 height ratio)
- Shared x-axis
- Residuals in barn units
- Zero reference line

### 4. Combined Visualization

Use all features together:

```python
ax_main, ax_res = result.plot_total_xs(
    split_phases=True,
    plot_residuals=True
)
```

## New CrossSection Methods

Added two utility methods to support the enhancements:

### `get_phase_xs(wl, phase)`

Get cross section for a specific phase:

```python
xs = nbragg.CrossSection(
    alpha="Fe_sg229_Iron-alpha.ncmat",
    gamma="Fe_sg225_Iron-gamma.ncmat"
)

# Get cross section for just the alpha phase
xs_alpha = xs.get_phase_xs(wavelength, 'alpha')
```

### `get_atomic_density()`

Get atomic number density:

```python
n = xs.get_atomic_density()  # atoms/(barn*cm)
```

## Usage Examples

### Basic Usage (Now Shows Data)

```python
import nbragg

# Create model and fit
xs = nbragg.CrossSection(iron="Fe_sg229_Iron-alpha.ncmat")
model = nbragg.TransmissionModel(xs, vary_basic=True)
result = model.fit(data, wlmin=1.0, wlmax=6.0)

# Plot - now shows both data and model!
result.plot_total_xs()
```

### Multi-Phase with Phase Contributions

```python
# Multi-phase material
xs = (nbragg.CrossSection(alpha="Fe_sg229_Iron-alpha.ncmat") * 0.7 +
      nbragg.CrossSection(gamma="Fe_sg225_Iron-gamma.ncmat") * 0.3)

# Fit with varying weights
model = nbragg.TransmissionModel(xs, vary_basic=True, vary_weights=True)
result = model.fit(data, wlmin=1.0, wlmax=6.0)

# Show individual phase contributions
result.plot_total_xs(split_phases=True)
```

### Complete Analysis Plot

```python
# Create comprehensive 2-panel plot
ax_main, ax_res = result.plot_total_xs(
    split_phases=True,        # Show all phases
    plot_residuals=True,      # Add residuals panel
    plot_bg=True,             # Show background
    title="Steel Texture Analysis"
)
```

### Custom Styling

```python
# Customize appearance
result.plot_total_xs(
    split_phases=True,
    plot_residuals=True,
    color="darkblue",         # Model line color
    title="My Custom Title"
)
```

## Technical Details

### Data Conversion

Transmission to cross section conversion:

```
T = norm * exp(-σ * thickness * n) * (1 - bg) + bg

Solving for σ:
σ = -ln((T - bg) / (norm * (1 - bg))) / (thickness * n)
```

Error propagation:
```
Δσ ≈ |σ/T| * ΔT
```

### Phase Cross Section Calculation

For each phase `i`:
```
xs_phase_i = get_phase_xs(wavelength, phase_i)
weighted_xs_i = xs_phase_i * weight_i
```

Total cross section:
```
xs_total = Σ weighted_xs_i
```

### Residuals Calculation

```
residuals = data_xs - xs_total
```

## Files Modified

### [models.py](src/nbragg/models.py)

**Lines 1522-1765:** Enhanced `plot_total_xs` method
- Added `split_phases` parameter
- Added `plot_residuals` parameter
- Data extraction and conversion
- Phase contribution plotting
- Residuals panel
- Updated docstring

**Lines 188:** Fixed SyntaxWarning (escaped backslashes)

**Lines 1622-1624:** Fixed error calculation (added `np.abs`)

### [cross_section.py](src/nbragg/cross_section.py)

**Lines 957-974:** Added `get_phase_xs` method

**Lines 976-999:** Added `get_atomic_density` method

## Backward Compatibility

All changes are backward compatible:
- Default behavior unchanged (no phases split, no residuals)
- New parameters are optional
- Existing code continues to work

## Visual Comparison

### Before
- Only model curve visible
- No indication of data quality
- No phase separation

### After
- ✅ Data points with error bars
- ✅ Individual phase contributions (optional)
- ✅ Chi-squared in legend
- ✅ Residuals panel (optional)
- ✅ Phase weight percentages

## Example Output

When using `split_phases=True` and `plot_residuals=True`, you get a plot similar to your example:

```python
result.plot_total_xs(split_phases=True, plot_residuals=True)
```

Creates:
- **Top panel:**
  - Data points (light purple with error bars)
  - Total cross section (black line)
  - Individual phases (colored lines from turbo colormap)
  - Phase labels with percentages (e.g., "α: 65.2%")
  - Chi-squared in legend

- **Bottom panel:**
  - Residuals (data - model)
  - Zero reference line
  - Same wavelength axis

## Testing

All features tested and working:
- ✅ Basic plot shows data
- ✅ Split phases displays all phases with colors
- ✅ Residuals panel works correctly
- ✅ Combined features work together
- ✅ Error bars calculated correctly
- ✅ Phase percentages displayed properly

The enhancement makes `plot_total_xs` a powerful tool for analyzing multi-phase Bragg edge data!

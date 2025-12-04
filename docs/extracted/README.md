# nBragg Crystal Orientation Documentation Package

This package contains comprehensive documentation for specifying crystal orientations in nBragg, ready for integration into your ReadTheDocs site.

## Files Included

### 1. `crystal_orientation_guide.rst` (Main Documentation)
**Purpose**: Comprehensive guide covering all aspects of crystal orientation specification

**Contents**:
- Overview of coordinate system conventions
- Mathematical framework (Euler angles, rotation matrices)
- Practical examples for different crystal systems
- Mosaicity and crystal imperfections
- Crystal rotation for Bragg dip alignment
- Visualization and interpretation guidelines
- Best practices and common pitfalls

**Target audience**: Users who need detailed understanding of the theory and implementation

---

### 2. `orientation_quick_reference.rst` (Quick Reference Card)
**Purpose**: Condensed reference for quick lookup during analysis

**Contents**:
- Coordinate system summary
- Common orientations table
- Key equations in compact form
- Troubleshooting checklist
- Code snippets for common tasks
- Mosaicity parameters reference

**Target audience**: Experienced users who need quick reminders

---

### 3. `orientation_examples.py` (Python Examples)
**Purpose**: Executable Python code demonstrating orientation specification

**Contents**:
- 7 complete working examples:
  1. Simple cubic orientation
  2. Rotated crystal (45°)
  3. Arbitrary orientation
  4. Rotation for alignment
  5. Orthogonality checking
  6. Hexagonal crystal system
  7. Mosaicity considerations
- Utility functions for rotation matrices
- Orthogonality verification tools

**Target audience**: All users, especially those who learn by example

---

## Integration into ReadTheDocs

### Recommended Documentation Structure

```
docs/
├── index.rst
├── user_guide/
│   ├── crystal_orientation_guide.rst      # Main guide
│   └── orientation_quick_reference.rst    # Quick ref
├── examples/
│   └── orientation_examples.py            # Code examples
└── api/
    └── ...
```

### Adding to Table of Contents

In your main `index.rst` or relevant section, add:

```rst
.. toctree::
   :maxdepth: 2
   :caption: User Guide:

   user_guide/crystal_orientation_guide
   user_guide/orientation_quick_reference

.. toctree::
   :maxdepth: 1
   :caption: Examples:

   examples/orientation_examples
```

### Sphinx Configuration

Ensure your `conf.py` includes:

```python
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',  # For LaTeX equations
    'sphinx.ext.viewcode',
]

# If using custom CSS for tables
html_static_path = ['_static']
```

---

## Key Features of This Documentation

### 1. **Progressive Complexity**
- Starts with simple concepts
- Builds to advanced topics
- Multiple entry points for different skill levels

### 2. **Mathematical Rigor**
- Proper LaTeX formatting for all equations
- Consistent notation throughout
- References to original papers (Matthies, NCrystal, etc.)

### 3. **Practical Focus**
- Real-world examples
- Working Python code
- Troubleshooting guidance
- Best practices from experience

### 4. **Visual Support**
- References to coordinate system figures
- Descriptions of Bragg circles visualization
- Table-based quick references

---

## Customization Notes

### Where Your Figure Fits

The documentation references a figure that should be included:

**Location in text**: Section "Visualization and Interpretation"

**Figure description needed**:
```rst
.. figure:: ../figures/nBraggCoords_cropped_1.pdf
   :width: 70%
   :align: center

   nBragg cartesian coordinate system showing incident neutron direction,
   crystal orientation vectors, rotation axes, and Bragg circles.
```

**Key elements your figure should show**:
- (x, y, z) laboratory coordinate system
- Incident beam k_i along +z
- [hkl]_z and [hkl]_y orientation vectors
- Rotation axes and angles (φ, θ)
- Bragg circles at different wavelengths
- xz-plane reference
- Mosaicity distribution around normal vector

### Modifying for Your API

The examples use generic function calls. Update these to match your actual nBragg API:

```python
# Generic example (in current docs):
# Create or update crystal with orientation

# Update to your actual API:
crystal = nbragg.Crystal(
    material='Fe',
    orientation_z=[0, 0, 1],
    orientation_y=[0, 1, 0],
    mosaicity_fwhm=0.5
)
```

### Adding Crystal System Specifics

Consider adding subsections for specific crystal systems users commonly work with:
- FCC metals (Cu, Al, Ni)
- BCC metals (Fe, W, Mo)
- HCP metals (Ti, Zr, Mg)
- Each with recommended orientations and examples

---

## Testing the Documentation

### 1. Build Locally
```bash
cd docs/
make html
# Check _build/html/crystal_orientation_guide.html
```

### 2. Verify Equations
- All LaTeX should render properly with MathJax
- Check that matrices display correctly
- Verify inline vs display math

### 3. Test Examples
```bash
python orientation_examples.py
```
Should run without errors and produce informative output

### 4. Check Links
Ensure cross-references work:
- `:doc:` links to other pages
- `:ref:` links to sections
- External citations

---

## Future Enhancements

Consider adding:

1. **Interactive Jupyter Notebook**
   - Visualize rotations in 3D
   - Interactive orientation selection
   - Real-time Bragg circle calculations

2. **Video Tutorials**
   - Screencast of orientation specification workflow
   - 3D visualization of coordinate transformations

3. **FAQ Section**
   - Common orientation questions from users
   - Troubleshooting specific crystal systems

4. **Advanced Topics**
   - Texture analysis workflows
   - Orientation distribution functions
   - Multi-phase materials

---

## Citations and References

The documentation references these key papers:

1. **Matthies (1987)**: Euler angle conventions
2. **Kittelmann et al. (2021)**: NCrystal implementation and mosaicity
3. **Schmitt et al. (2023)**: NCrystal orientation handling
4. **Cai et al. (2020)**: Bragg-edge imaging theory

**Action needed**: Update these citations in your bibliography:
```rst
.. [Matthies1987] Matthies, S. (1987). On the diffracting volumes...
.. [Kittelmann2021] Kittelmann, T., et al. (2021). ...
.. [Schmitt2023] Schmitt, B., et al. (2023). ...
.. [Cai2020] Cai, W., et al. (2020). ...
```

---

## Support and Feedback

This documentation should address most user questions about crystal orientation. Consider:

1. Adding a "Was this helpful?" feedback widget
2. Tracking which sections users visit most
3. Collecting user questions to improve docs
4. Creating a discussion forum for orientation-related questions

---

## License and Attribution

Include appropriate licensing:
```rst
License
-------
This documentation is part of the nBragg package.

[Your license terms here]

Contributors
------------
- [Your name] - Original implementation
- [Colleague's name] - Coordinate system figure and guidance
- nBragg community - Feedback and suggestions
```

---

## Quick Start for Documentation Authors

1. **Copy files to your docs directory**
2. **Add figure reference (nBraggCoords_cropped_1.pdf)**
3. **Update API examples to match your code**
4. **Add to table of contents**
5. **Build and test**
6. **Push to ReadTheDocs**

---

## Contact

For questions about this documentation package:
- Check the examples first
- Refer to the quick reference card
- Consult the main guide for theory
- File an issue if something is unclear

---

**Version**: 1.0
**Date**: December 2024
**Status**: Ready for integration

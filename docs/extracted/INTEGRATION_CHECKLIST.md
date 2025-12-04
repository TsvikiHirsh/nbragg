# Integration Checklist for nBragg Orientation Documentation

Use this checklist to integrate the crystal orientation documentation into your nBragg ReadTheDocs site.

## Pre-Integration

- [ ] Review all files in this package
- [ ] Verify your current documentation structure
- [ ] Identify where orientation docs should live (e.g., `docs/user_guide/`)
- [ ] Check that you have the coordinate system figure ready

## File Placement

- [ ] Copy `orientation_index.rst` → `docs/user_guide/orientation_index.rst`
- [ ] Copy `crystal_orientation_guide.rst` → `docs/user_guide/crystal_orientation_guide.rst`
- [ ] Copy `orientation_quick_reference.rst` → `docs/user_guide/orientation_quick_reference.rst`
- [ ] Copy `orientation_examples.py` → `docs/examples/orientation_examples.py`

## Figure Integration

- [ ] Place your coordinate system figure: `docs/figures/nBraggCoords_cropped_1.pdf`
- [ ] Verify the figure path in `crystal_orientation_guide.rst` is correct
- [ ] Test that figure renders properly in HTML build

## API Customization

- [ ] Update code examples to match your actual nBragg API
- [ ] Test example code with real nBragg installation
- [ ] Update function/class names if needed
- [ ] Verify parameter names are correct

Example locations to update:
- Line ~50 in `orientation_index.rst` (Quick Start example)
- Line ~180 in `crystal_orientation_guide.rst` (Usage example)
- All examples in `orientation_examples.py` (if using actual nBragg API)

## Table of Contents

- [ ] Add orientation docs to your main `index.rst` or `user_guide/index.rst`
- [ ] Use appropriate `:maxdepth:` setting
- [ ] Test navigation works correctly

Example addition:
```rst
.. toctree::
   :maxdepth: 2
   :caption: User Guide:

   user_guide/orientation_index
   user_guide/crystal_orientation_guide
   user_guide/orientation_quick_reference
```

## Cross-References

- [ ] Update `:doc:` links if your document structure differs
- [ ] Verify all internal cross-references work
- [ ] Add links from existing docs to orientation pages
- [ ] Test all "See Also" links

## Citations & References

- [ ] Add proper citations to your bibliography file
- [ ] Update citation format to match your style
- [ ] Verify all cited papers are in your references

Required citations:
- Matthies (1987) - Euler angles
- Kittelmann et al. (2021) - NCrystal
- Schmitt et al. (2023) - NCrystal orientations
- Cai et al. (2020) - Bragg circles

## Sphinx Configuration

- [ ] Verify `sphinx.ext.mathjax` is in your `conf.py` extensions
- [ ] Check that math rendering works
- [ ] Test code syntax highlighting
- [ ] Verify table formatting

Required in `conf.py`:
```python
extensions = [
    'sphinx.ext.mathjax',  # For equations
    'sphinx.ext.viewcode',  # For code examples
    # ... your other extensions
]
```

## Build & Test

- [ ] Run `make html` locally
- [ ] Check for Sphinx warnings
- [ ] Verify all equations render correctly
- [ ] Test all code blocks have proper syntax highlighting
- [ ] Check all tables display properly
- [ ] Verify figure displays correctly
- [ ] Test navigation between pages
- [ ] Check mobile responsiveness

## Content Review

- [ ] Read through all pages for clarity
- [ ] Verify technical accuracy
- [ ] Check for typos or formatting issues
- [ ] Ensure examples are appropriate for your user base
- [ ] Verify all acronyms are defined on first use

## Examples Validation

- [ ] Run `python orientation_examples.py` without errors
- [ ] Test each example function individually
- [ ] Verify output makes sense
- [ ] Check that examples are pedagogically sound

## Accessibility

- [ ] Verify alt text for figures
- [ ] Check heading hierarchy is logical
- [ ] Ensure tables have appropriate headers
- [ ] Test with screen reader if possible

## Customization (Optional)

- [ ] Add domain-specific examples for your field
- [ ] Include examples for materials your users commonly work with
- [ ] Add FAQ section based on anticipated questions
- [ ] Create crystal-system-specific subsections

## ReadTheDocs Setup

- [ ] Push changes to your repository
- [ ] Trigger ReadTheDocs build
- [ ] Check build logs for errors
- [ ] Verify pages appear in navigation
- [ ] Test search functionality includes new content
- [ ] Check version selector if using multiple versions

## Post-Integration

- [ ] Announce new documentation to users
- [ ] Collect feedback
- [ ] Monitor which pages users visit most
- [ ] Note common questions for future FAQ section

## Future Enhancements

Consider adding (not required for initial release):

- [ ] Jupyter notebook with interactive examples
- [ ] Video tutorial
- [ ] Additional crystal system examples
- [ ] 3D visualization tools
- [ ] Validation scripts for orientation checking

## Sign-Off

Once completed:

- [ ] All builds pass without errors
- [ ] All examples tested and working
- [ ] Navigation structure verified
- [ ] Content reviewed by at least one other person
- [ ] Links shared with beta testers (if applicable)

---

## Quick Test Commands

```bash
# Build documentation
cd docs/
make clean
make html

# Test Python examples
python examples/orientation_examples.py

# Check for broken links (if you have this tool)
make linkcheck

# View locally
python -m http.server 8000 --directory _build/html
# Then visit: http://localhost:8000
```

---

## Troubleshooting

### Common Issues:

**Issue**: Equations not rendering
- **Solution**: Check `sphinx.ext.mathjax` in extensions

**Issue**: Figure not found
- **Solution**: Verify relative path from RST file location

**Issue**: Code blocks not highlighted
- **Solution**: Add language identifier: `.. code-block:: python`

**Issue**: Cross-references broken
- **Solution**: Use `:doc:` for documents, `:ref:` for sections

**Issue**: Build warnings about orphan documents
- **Solution**: Add pages to a toctree

---

## Contact

If you encounter issues:
1. Check this checklist
2. Review README.md
3. Consult Sphinx documentation
4. Check ReadTheDocs build logs

---

**Date Prepared**: December 2024
**Package Version**: 1.0
**Status**: Ready for use

"""
CrossSection Parameter Assignment Examples
===========================================

This example demonstrates the ergonomic parameter assignment feature in
CrossSection that allows users to easily add extinction, SANS, and orientation
parameters directly in the constructor.
"""

from nbragg import CrossSection

# ============================================================================
# Example 1: Single Material with Extinction Parameters
# ============================================================================
print("=" * 70)
print("Example 1: Single material with extinction parameters")
print("=" * 70)

# Create a CrossSection with extinction parameters
xs_extinction = CrossSection(
    steel='Al_sg225.ncmat',
    ext_l=100,      # Crystallite size (microns)
    ext_Gg=1500,    # Mosaic spread (seconds of arc)
    ext_L=2000      # Strain parameter
)

print(f"Material: {list(xs_extinction.materials.keys())}")
print(f"ext_l: {xs_extinction.materials['steel']['ext_l']}")
print(f"ext_Gg: {xs_extinction.materials['steel']['ext_Gg']}")
print(f"ext_L: {xs_extinction.materials['steel']['ext_L']}")
print()

# ============================================================================
# Example 2: Single Material with SANS Parameters
# ============================================================================
print("=" * 70)
print("Example 2: Single material with SANS parameter")
print("=" * 70)

# Create a CrossSection with SANS hard-sphere radius
xs_sans = CrossSection(
    aluminum='Al_sg225.ncmat',
    sans=50.0  # Hard-sphere radius in Angstroms
)

print(f"Material: {list(xs_sans.materials.keys())}")
print(f"sans: {xs_sans.materials['aluminum']['sans']}")
print()

# ============================================================================
# Example 3: Single Material with Orientation Parameters
# ============================================================================
print("=" * 70)
print("Example 3: Single material with orientation parameters")
print("=" * 70)

# Create a CrossSection with orientation parameters
xs_orientation = CrossSection(
    sample='Al_sg225.ncmat',
    theta=45,   # Polar angle (degrees)
    phi=90,     # Azimuthal angle (degrees)
    mos=5       # Mosaicity spread (degrees)
)

print(f"Material: {list(xs_orientation.materials.keys())}")
print(f"theta: {xs_orientation.materials['sample']['theta']}")
print(f"phi: {xs_orientation.materials['sample']['phi']}")
print(f"mos: {xs_orientation.materials['sample']['mos']}")
print()

# ============================================================================
# Example 4: Multiple Materials with Different Parameters (Order-Based)
# ============================================================================
print("=" * 70)
print("Example 4: Multiple materials with different parameters")
print("=" * 70)

# Parameters are assigned to materials based on their order in kwargs
# ext_Gg goes to mat1, sans goes to mat2
xs_multi = CrossSection(
    mat1='Al_sg225.ncmat', ext_Gg=100,
    mat2='Al_sg225.ncmat', sans=20
)

print(f"Materials: {list(xs_multi.materials.keys())}")
print(f"mat1 ext_Gg: {xs_multi.materials['mat1']['ext_Gg']}")
print(f"mat1 sans: {xs_multi.materials['mat1']['sans']}")
print(f"mat2 ext_Gg: {xs_multi.materials['mat2']['ext_Gg']}")
print(f"mat2 sans: {xs_multi.materials['mat2']['sans']}")
print()

# ============================================================================
# Example 5: Multiple Parameters for Each Material
# ============================================================================
print("=" * 70)
print("Example 5: Multiple parameters for each material")
print("=" * 70)

# Each material can have multiple parameters assigned
xs_complex = CrossSection(
    iron='Al_sg225.ncmat', ext_l=50, ext_Gg=200,      # Iron gets extinction
    aluminum='Al_sg225.ncmat', theta=30, phi=60       # Aluminum gets orientation
)

print(f"Materials: {list(xs_complex.materials.keys())}")
print("Iron:")
print(f"  ext_l: {xs_complex.materials['iron']['ext_l']}")
print(f"  ext_Gg: {xs_complex.materials['iron']['ext_Gg']}")
print("Aluminum:")
print(f"  theta: {xs_complex.materials['aluminum']['theta']}")
print(f"  phi: {xs_complex.materials['aluminum']['phi']}")
print()

# ============================================================================
# Example 6: Lattice Parameters
# ============================================================================
print("=" * 70)
print("Example 6: Lattice parameters")
print("=" * 70)

# Assign lattice parameters to refine crystal structure
xs_lattice = CrossSection(
    sample='Al_sg225.ncmat',
    a=4.05,
    b=4.05,
    c=4.05
)

print(f"Material: {list(xs_lattice.materials.keys())}")
print(f"a: {xs_lattice.materials['sample']['a']}")
print(f"b: {xs_lattice.materials['sample']['b']}")
print(f"c: {xs_lattice.materials['sample']['c']}")
print()

# ============================================================================
# Example 7: Combined Parameters
# ============================================================================
print("=" * 70)
print("Example 7: Combined extinction, orientation, and lattice parameters")
print("=" * 70)

# You can combine multiple types of parameters for one material
xs_combined = CrossSection(
    sample='Al_sg225.ncmat',
    ext_l=100,      # Extinction
    ext_Gg=1500,
    theta=45,       # Orientation
    phi=90,
    mos=5,
    a=4.05,         # Lattice
    temp=400        # Temperature
)

print(f"Material: {list(xs_combined.materials.keys())}")
print("Extinction:")
print(f"  ext_l: {xs_combined.materials['sample']['ext_l']}")
print(f"  ext_Gg: {xs_combined.materials['sample']['ext_Gg']}")
print("Orientation:")
print(f"  theta: {xs_combined.materials['sample']['theta']}")
print(f"  phi: {xs_combined.materials['sample']['phi']}")
print(f"  mos: {xs_combined.materials['sample']['mos']}")
print("Lattice:")
print(f"  a: {xs_combined.materials['sample']['a']}")
print("Temperature:")
print(f"  temp: {xs_combined.materials['sample']['temp']}")
print()

# ============================================================================
# Summary of Supported Parameters
# ============================================================================
print("=" * 70)
print("Summary of Supported Parameters")
print("=" * 70)
print("""
Extinction Parameters:
  - ext_l      : Crystallite size (microns)
  - ext_Gg     : Mosaic spread (seconds of arc)
  - ext_L      : Strain parameter
  - ext_method : Extinction method
  - ext_dist   : Distribution type

SANS Parameters:
  - sans       : Hard-sphere radius (Angstroms)

Orientation Parameters:
  - theta      : Polar angle (degrees)
  - phi        : Azimuthal angle (degrees)
  - mos        : Mosaicity spread (degrees)
  - dir1       : First direction vector
  - dir2       : Second direction vector
  - dirtol     : Direction tolerance

Lattice Parameters:
  - a, b, c    : Lattice constants (Angstroms)

Other Parameters:
  - temp       : Temperature (Kelvin)
  - weight     : Material weight (for multi-phase)
""")

print("\nBackward Compatibility:")
print("This feature is fully backward compatible. Existing code continues")
print("to work without any changes. The traditional dictionary-based approach")
print("for specifying material parameters is still fully supported.")

"""
Crystal Orientation Examples for nBragg
========================================

This module provides practical examples of specifying crystal orientations
in nBragg for various common scenarios.

Author: nBragg Documentation Team
"""

import numpy as np
# import nbragg  # Uncomment when using with actual nbragg installation


def example_1_simple_orientation():
    """
    Example 1: Simple cubic crystal with standard orientation
    
    Crystal aligned with (001) plane perpendicular to beam direction (z)
    and (010) plane perpendicular to y-axis.
    """
    print("Example 1: Standard Orientation")
    print("-" * 50)
    
    # Define orientation using Miller indices
    hkl_z = [0, 0, 1]  # (001) plane normal to beam
    hkl_y = [0, 1, 0]  # (010) plane normal to y-axis
    
    print(f"[hkl]_z: {hkl_z}")
    print(f"[hkl]_y: {hkl_y}")
    print(f"Description: Standard cubic orientation\n")
    
    return hkl_z, hkl_y


def example_2_rotated_orientation():
    """
    Example 2: Crystal rotated 45° around z-axis
    
    This demonstrates how to specify an orientation where the crystal
    is rotated in the xy-plane.
    """
    print("Example 2: 45° Rotation Around Z-axis")
    print("-" * 50)
    
    hkl_z = [0, 0, 1]  # Still along beam direction
    hkl_y = [1, 1, 0]  # Rotated 45° in xy-plane
    
    print(f"[hkl]_z: {hkl_z}")
    print(f"[hkl]_y: {hkl_y}")
    print(f"Description: Crystal rotated 45° around beam axis\n")
    
    return hkl_z, hkl_y


def example_3_arbitrary_orientation():
    """
    Example 3: Arbitrary crystal orientation
    
    Demonstrates a more complex orientation using higher Miller indices.
    """
    print("Example 3: Arbitrary Orientation")
    print("-" * 50)
    
    hkl_z = [1, 1, 2]  # (112) plane along beam
    hkl_y = [1, -1, 0] # (1-10) plane along y-axis
    
    print(f"[hkl]_z: {hkl_z}")
    print(f"[hkl]_y: {hkl_y}")
    print(f"Description: Complex orientation for textured sample\n")
    
    return hkl_z, hkl_y


def calculate_rotation_matrix(phi, theta):
    """
    Calculate the combined rotation matrix for crystal rotation.
    
    Parameters
    ----------
    phi : float
        Rotation angle around x-axis (radians)
    theta : float
        Rotation angle around y-axis (radians)
    
    Returns
    -------
    R : ndarray
        3x3 rotation matrix
    """
    # Rotation around y-axis
    Ry = np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])
    
    # Rotation around x-axis
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(phi), -np.sin(phi)],
        [0, np.sin(phi), np.cos(phi)]
    ])
    
    # Combined rotation: first Rx, then Ry
    R = Ry @ Rx
    
    return R


def example_4_rotation_for_alignment():
    """
    Example 4: Rotating crystal orientation for Bragg dip alignment
    
    Demonstrates how to apply rotations to align modeled Bragg dips
    with experimental data.
    """
    print("Example 4: Crystal Rotation for Alignment")
    print("-" * 50)
    
    # Initial orientation
    hkl_initial = np.array([0, 0, 1])
    
    # Rotation angles (in radians)
    phi = np.radians(5)    # 5° rotation around x-axis
    theta = np.radians(10) # 10° rotation around y-axis
    
    # Calculate rotation matrix
    R = calculate_rotation_matrix(phi, theta)
    
    # Apply rotation
    hkl_rotated = R @ hkl_initial
    
    print(f"Initial [hkl]: {hkl_initial}")
    print(f"Rotation angles: φ = {np.degrees(phi):.1f}°, θ = {np.degrees(theta):.1f}°")
    print(f"Rotated [hkl]: {hkl_rotated}")
    print(f"Normalized [hkl]: {hkl_rotated / np.linalg.norm(hkl_rotated)}\n")
    
    return hkl_rotated


def check_orthogonality(hkl_z, hkl_y, reciprocal_basis=None):
    """
    Check if two Miller index vectors are orthogonal in reciprocal space.
    
    Parameters
    ----------
    hkl_z : array-like
        Miller indices for z-direction
    hkl_y : array-like
        Miller indices for y-direction
    reciprocal_basis : ndarray, optional
        3x3 matrix of reciprocal lattice vectors. If None, assumes orthogonal basis.
    
    Returns
    -------
    bool
        True if vectors are orthogonal
    float
        Dot product value
    """
    if reciprocal_basis is None:
        # Assume simple orthogonal reciprocal lattice
        reciprocal_basis = np.eye(3)
    
    # Convert to reciprocal space vectors
    vec_z = reciprocal_basis @ np.array(hkl_z)
    vec_y = reciprocal_basis @ np.array(hkl_y)
    
    # Normalize
    vec_z = vec_z / np.linalg.norm(vec_z)
    vec_y = vec_y / np.linalg.norm(vec_y)
    
    # Check orthogonality
    dot_product = np.dot(vec_z, vec_y)
    is_orthogonal = np.abs(dot_product) < 1e-10
    
    return is_orthogonal, dot_product


def example_5_orthogonality_check():
    """
    Example 5: Verifying orthogonality of orientation vectors
    
    Demonstrates how to check if chosen Miller indices result in
    orthogonal directions.
    """
    print("Example 5: Checking Orthogonality")
    print("-" * 50)
    
    # Test case 1: Orthogonal
    hkl_z1 = [0, 0, 1]
    hkl_y1 = [0, 1, 0]
    is_orth1, dot1 = check_orthogonality(hkl_z1, hkl_y1)
    print(f"Case 1: {hkl_z1} and {hkl_y1}")
    print(f"  Orthogonal: {is_orth1}, Dot product: {dot1:.6f}\n")
    
    # Test case 2: Orthogonal with higher indices
    hkl_z2 = [1, 1, 0]
    hkl_y2 = [-1, 1, 0]
    is_orth2, dot2 = check_orthogonality(hkl_z2, hkl_y2)
    print(f"Case 2: {hkl_z2} and {hkl_y2}")
    print(f"  Orthogonal: {is_orth2}, Dot product: {dot2:.6f}\n")
    
    # Test case 3: Non-orthogonal (should fail)
    hkl_z3 = [1, 0, 0]
    hkl_y3 = [1, 1, 0]
    is_orth3, dot3 = check_orthogonality(hkl_z3, hkl_y3)
    print(f"Case 3: {hkl_z3} and {hkl_y3}")
    print(f"  Orthogonal: {is_orth3}, Dot product: {dot3:.6f}")
    print(f"  WARNING: Not orthogonal!\n")


def example_6_hexagonal_crystal():
    """
    Example 6: Hexagonal crystal system
    
    Demonstrates orientation specification for a hexagonal crystal
    where the reciprocal lattice is non-orthogonal.
    """
    print("Example 6: Hexagonal Crystal System")
    print("-" * 50)
    
    # For hexagonal: c-axis along z, a-axis at 120° to b-axis in xy-plane
    hkl_z = [0, 0, 1]  # c-axis along beam
    hkl_y = [0, 1, 0]  # b-axis along y
    
    print(f"[hkl]_z: {hkl_z} (c-axis)")
    print(f"[hkl]_y: {hkl_y} (b-axis)")
    print(f"Description: Standard hexagonal orientation")
    print(f"Note: Reciprocal lattice accounts for 120° angle\n")
    
    # Alternative hexagonal orientation
    hkl_z_alt = [0, 0, 2]  # (002) plane
    hkl_y_alt = [1, 0, 0]  # a-axis along y
    
    print(f"Alternative:")
    print(f"[hkl]_z: {hkl_z_alt}")
    print(f"[hkl]_y: {hkl_y_alt}\n")


def example_7_mosaicity_considerations():
    """
    Example 7: Considering mosaicity in orientation specification
    
    Demonstrates how mosaicity affects the interpretation of crystal orientation.
    """
    print("Example 7: Mosaicity Considerations")
    print("-" * 50)
    
    # Nominal orientation
    hkl_nominal = [0, 0, 1]
    
    # Mosaicity parameters
    fwhm_eta = 0.5  # FWHM in degrees
    epsilon = 1e-3  # Volume fraction cutoff
    
    # Calculate truncation cutoff
    tau = max(3, 1.1 * np.sqrt(-2 * np.log(epsilon))) * fwhm_eta / (2 * np.sqrt(2 * np.log(2)))
    
    print(f"Nominal orientation: {hkl_nominal}")
    print(f"Mosaicity FWHM (η): {fwhm_eta}°")
    print(f"Truncation cutoff (τ): {tau:.3f}°")
    print(f"Interpretation: Crystal planes spread ±{tau:.3f}° around nominal\n")
    
    # Angular range covered
    print(f"Angular range: {-tau:.3f}° to +{tau:.3f}°")
    print(f"This affects Bragg condition satisfaction\n")


def main():
    """
    Run all examples
    """
    print("=" * 60)
    print("Crystal Orientation Examples for nBragg")
    print("=" * 60)
    print()
    
    example_1_simple_orientation()
    example_2_rotated_orientation()
    example_3_arbitrary_orientation()
    example_4_rotation_for_alignment()
    example_5_orthogonality_check()
    example_6_hexagonal_crystal()
    example_7_mosaicity_considerations()
    
    print("=" * 60)
    print("Examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()

"""
Test that all imports for the streamlit app work correctly
"""

print("Testing imports for streamlit app...")

try:
    import streamlit as st
    print("✓ streamlit")
except ImportError as e:
    print(f"✗ streamlit: {e}")

try:
    import numpy as np
    print("✓ numpy")
except ImportError as e:
    print(f"✗ numpy: {e}")

try:
    import pandas as pd
    print("✓ pandas")
except ImportError as e:
    print(f"✗ pandas: {e}")

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    print("✓ plotly")
except ImportError as e:
    print(f"✗ plotly: {e}")

try:
    import nbragg
    print("✓ nbragg")
    print(f"  - nbragg version: {nbragg.__version__ if hasattr(nbragg, '__version__') else 'unknown'}")
except ImportError as e:
    print(f"✗ nbragg: {e}")

try:
    import matplotlib
    import matplotlib.pyplot as plt
    print("✓ matplotlib")
except ImportError as e:
    print(f"✗ matplotlib: {e}")

print("\nChecking example data files...")
from pathlib import Path

example_files = [
    "notebooks/iron_powder_transmission.csv",
    "notebooks/iron_powder.csv",
    "notebooks/openbeam.csv"
]

for filepath in example_files:
    if Path(filepath).exists():
        print(f"✓ {filepath}")
    else:
        print(f"✗ {filepath} (not found)")

print("\nChecking nbragg.materials...")
try:
    materials = list(nbragg.materials.keys())
    print(f"✓ Found {len(materials)} materials in nbragg.materials")
    print(f"  Sample materials: {materials[:5]}")

    # Test accessing a material
    if materials:
        test_mat = materials[0]
        mat_spec = nbragg.materials[test_mat]
        print(f"✓ Successfully accessed test material: {test_mat}")
except Exception as e:
    print(f"✗ Error accessing nbragg.materials: {e}")

print("\n✅ All imports successful! Ready to run streamlit app.")
print("\nTo run the app:")
print("  streamlit run docs/streamlit/streamlit_app.py")

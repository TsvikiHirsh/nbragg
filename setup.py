from setuptools import setup, Extension, find_packages
import pybind11

# Define the extension module
ext_modules = [
    Extension(
        "nres._integrate_xs",
        sources=["src/bindings.cpp", "src/integrate_xs.cpp"],  # Source files
        include_dirs=[
            pybind11.get_include(),  # Include Pybind11 headers
            "include"  # Add the 'include' folder where your header files are located
        ],
    )
]

# Use setup to build the extension module
setup(
    name="nres",
    version="0.1.0",
    author="Tsviki Y. Hirsh",
    author_email="tsviki.hirsh@gmail.com",
    description="Simple yet powerful package for neutron resonance fitting",
    packages=find_packages(where="src"),  # Look for packages in src/
    package_dir={"": "src"},  # Tell setuptools that packages are in src/
    include_package_data=True,  # To include non-python files if any
    python_requires=">=3.8",
    ext_modules=ext_modules,
    zip_safe=False,  # Avoids problems with C++ extensions and zip imports
)
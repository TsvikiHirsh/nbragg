============
Installation
============

nbragg requires Python 3.9 or later.

From PyPI
---------

Install the latest stable version directly from PyPI:

.. code-block:: bash

    pip install nbragg

From Source
-----------

To install from the source repository:

.. code-block:: bash

    git clone https://github.com/TsvikiHirsh/nbragg
    cd nbragg
    pip install .

Dependencies
------------

nbragg requires the following packages, which are installed automatically by pip:

- NCrystal (cross-section calculations)
- lmfit (nonlinear fitting)
- numpy, scipy, pandas (numerics and data handling)
- matplotlib (plotting)
- spglib (crystal structure handling for SANS materials)
- joblib, tqdm (parallel fitting and progress bars)

Optional Dependencies
---------------------

**Extinction effects** require the CrysExtn NCrystal plugin:

.. code-block:: bash

    pip install git+https://github.com/XuShuqi7/ncplugin-CrysExtn

**Interactive widgets** (``model.interactive_plot()`` in Jupyter):

.. code-block:: bash

    pip install "nbragg[interactive]"

**Development tools** (running tests, pre-commit hooks):

.. code-block:: bash

    pip install "nbragg[dev]"

# nbragg tutorial notebooks

Three Jupyter notebooks cover the main nbragg workflows. Recommended order:

| # | Notebook | What it covers |
|---|----------|----------------|
| 1 | [nbragg_tutorial.ipynb](nbragg_tutorial.ipynb) | Getting started: loading data, defining cross-sections, transmission models, fitting, staged fits, SANS, extinction, oriented materials, and MTEX integration. |
| 2 | [Rietveld_in_nbragg_tutorial.ipynb](Rietveld_in_nbragg_tutorial.ipynb) | Rietveld-type analysis: staged, parametric refinement of Bragg edge data, stage summaries, and convergence diagnostics. |
| 3 | [grouped_fits_tutorial.ipynb](grouped_fits_tutorial.ipynb) | Grouped/gridded data: loading spatially-resolved or multi-sample data, parallel fitting, parameter maps, and saving/loading grouped results. |

## Running the notebooks

The notebooks use the data files (`*.csv`, `*.ncmat`) in this directory, so run
them from here:

```bash
pip install nbragg jupyter
git clone https://github.com/TsvikiHirsh/nbragg
cd nbragg/notebooks
jupyter lab
```

The first and second tutorials use measured iron-powder and large-grain steel
transmission data included in this directory; the grouped-fitting tutorial
generates its own synthetic example data.

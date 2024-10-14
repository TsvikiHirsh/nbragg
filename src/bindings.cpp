#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "integrate_xs.h"

namespace py = pybind11;

PYBIND11_MODULE(_integrate_xs, m) {
    m.def("integrate_cross_section", &integrate_cross_section, 
          "Integrates cross-section data with optional response and flight path length.",
          py::arg("xs_energies"), 
          py::arg("xs_values"), 
          py::arg("energy_grid"), 
          py::arg("response") = std::vector<double>{0.,1.,0.} // Default response
    );
}

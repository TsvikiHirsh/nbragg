#ifndef INTEGRATE_XS_H
#define INTEGRATE_XS_H

#include <vector>
#include <string>  // For std::string

// Declare the function with default parameters
std::vector<double> integrate_cross_section(
    const std::vector<double>& xs_energies,  // Cross-section energies
    const std::vector<double>& xs_values,    // Cross-section values
    const std::vector<double>& energy_grid,  // User energy grid
    const std::vector<double>& response = std::vector<double>{0.,1.,0.}  // Default response vector
);

#endif // INTEGRATE_XS_H

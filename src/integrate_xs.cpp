#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

// Linear interpolation function
double linear_interp(const std::vector<double>& xs_energies, const std::vector<double>& xs_values, double energy) {
    size_t n = xs_energies.size();
    
    // If the energy is out of bounds, return 0 (or handle as needed)
    if (energy <= xs_energies.front()) return xs_values.front();
    if (energy >= xs_energies.back()) return xs_values.back();

    // Find the two points in the grid between which the energy lies
    for (size_t i = 0; i < n - 1; ++i) {
        if (energy >= xs_energies[i] && energy <= xs_energies[i + 1]) {
            double slope = (xs_values[i + 1] - xs_values[i]) / (xs_energies[i + 1] - xs_energies[i]);
            return xs_values[i] + slope * (energy - xs_energies[i]);
        }
    }
    
    return 0.0; // Default return in case of error
}

// Convolve the values with the kernel and return the result with the same size as the input
std::vector<double> convolve_with_kernel(const std::vector<double>& values, const std::vector<double>& kernel) {
    int values_size = values.size();
    int kernel_size = kernel.size();
    int pad_size = kernel_size / 2; // Half of the kernel size for symmetric padding

    // Create a padded version of the input values
    std::vector<double> padded_values(values_size + 2 * pad_size, 0.0);
    
    // Copy original values into the padded vector
    std::copy(values.begin(), values.end(), padded_values.begin() + pad_size);

    // Result vector will have the same size as the original values
    std::vector<double> result(values_size, 0.0);

    // Perform convolution
    for (int i = 0; i < values_size; ++i) {
        double conv_sum = 0.0;
        for (int j = 0; j < kernel_size; ++j) {
            conv_sum += padded_values[i + j] * kernel[j];
        }
        result[i] = conv_sum;
    }

    return result;
}

std::vector<double> integrate_cross_section(
    const std::vector<double>& xs_energies,
    const std::vector<double>& xs_values,
    const std::vector<double>& energy_grid,
    const std::vector<double>& kernel)
{
    // Handle the case when xs_values has only one element
    if (xs_values.size() == 1) {
        return std::vector<double>(energy_grid.size(), xs_values[0]);
    }

    // Calculate the number of bins to add based on the kernel length
    int num_bins_to_add = kernel.empty() ? 0 : kernel.size() - 1;

    // Extend the energy grid
    std::vector<double> extended_energy_grid = energy_grid;
    for (int i = 0; i < num_bins_to_add; ++i) {
        // Add prefix bin
        double new_prefix = extended_energy_grid.front() * std::pow(extended_energy_grid.front() / extended_energy_grid[1], 1);
        extended_energy_grid.insert(extended_energy_grid.begin(), new_prefix);
        
        // Add suffix bin
        double new_suffix = extended_energy_grid.back() * std::pow(extended_energy_grid.back() / extended_energy_grid[extended_energy_grid.size() - 2], 1);
        extended_energy_grid.push_back(new_suffix);
    }

    std::vector<double> integrated_values;

    // Perform integration over the extended energy grid
    for (size_t i = 0; i < extended_energy_grid.size() - 1; ++i) {
        double emin = extended_energy_grid[i];
        double emax = extended_energy_grid[i + 1];

        if (emin <= 0 || emax <= 0) {
            integrated_values.push_back(0.0);
            continue;
        }

        // Perform trapezoidal integration
        double integral = 0.5 * (linear_interp(xs_energies, xs_values, emin) + 
                                 linear_interp(xs_energies, xs_values, emax)) * (emax - emin);

        // Calculate the average cross-section for this bin
        double avg_xs = integral / (emax - emin);
        integrated_values.push_back(avg_xs);
    }

    // If a kernel is provided and its size is greater than 1, convolve the result
    if (kernel.size() > 1) {
        integrated_values = convolve_with_kernel(integrated_values, kernel);
    }

    // Trim the result to match the original energy grid size
    size_t expected_size = energy_grid.size();
    if (integrated_values.size() > expected_size) {
        size_t extra = integrated_values.size() - expected_size;
        size_t to_remove_start = extra / 2;
        size_t to_remove_end = extra - to_remove_start;
        integrated_values.erase(integrated_values.begin(), integrated_values.begin() + to_remove_start);
        integrated_values.erase(integrated_values.end() - to_remove_end, integrated_values.end());
    }

    // Ensure the result matches the original energy grid size
    if (integrated_values.size() != expected_size) {
        integrated_values.resize(expected_size, integrated_values.back());
    }

    // Replace NaN values with zero
    std::replace_if(integrated_values.begin(), integrated_values.end(), 
                    [](double val) { return std::isnan(val); }, 0.0);

    return integrated_values;
}
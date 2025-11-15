#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <functional>
#include <stdexcept>
#include "interpolate.h"


inline std::vector<double> adjust_bout(std::vector<double> inarray, int fs = 10) {
    // if data is available for 70% of the last second
    int remainder = inarray.size() % fs;
    if (remainder >= static_cast<int>(0.7 * fs)) {
        for (int i = 0; i < fs - remainder; ++i) {
            inarray.push_back(inarray.back());
        }
    }
    // otherwise, trim the data to the full second
    else {
       inarray.resize((inarray.size() / fs) * fs);
    }  return inarray;
}


std::pair<std::vector<double>, std::vector<double>> preprocess_bout(
    const std::vector<double>& t_bout,
    const std::vector<double>& x_bout,
    const std::vector<double>& y_bout,
    const std::vector<double>& z_bout,
    int fs = 10
) {
    if (t_bout.size() < 2 || x_bout.size() < 2 || y_bout.size() < 2 || z_bout.size() < 2) {
        return {{}, {}};
    }

    // Create resampled time vector
    double t_start = t_bout.front();
    double t_end = t_bout.back();
    std::vector<double> t_bout_interp;

    for (double t = 0; t <= (t_end - t_start); t += 1.0 / fs) {
        t_bout_interp.push_back(t + t_start);
    }

    // Interpolate x, y, z
    std::vector<double> x_bout_interp = interp1d(t_bout, x_bout, t_bout_interp);
    std::vector<double> y_bout_interp = interp1d(t_bout, y_bout, t_bout_interp);
    std::vector<double> z_bout_interp = interp1d(t_bout, z_bout, t_bout_interp);

    // Adjust bouts
    x_bout_interp = adjust_bout(x_bout_interp);
    y_bout_interp = adjust_bout(y_bout_interp);
    z_bout_interp = adjust_bout(z_bout_interp);

    // Number of full seconds of measurements
    size_t num_seconds = static_cast<size_t>(std::floor(x_bout_interp.size() / static_cast<double>(fs)));

    // Trim and decimate t_bout_interp
    size_t trim_length = num_seconds * fs;
    if (trim_length > x_bout_interp.size()) {
        trim_length = x_bout_interp.size();
    }
    t_bout_interp.resize(trim_length);
    std::vector<double> t_bout_decimated;
    for (size_t i = 0; i < t_bout_interp.size(); i += fs) {
        t_bout_decimated.push_back(t_bout_interp[i]);
    }
    t_bout_interp = std::move(t_bout_decimated);

    // Calculate vector magnitude vm_bout_interp
    std::vector<double> vm_bout_interp;
    vm_bout_interp.reserve(x_bout_interp.size());
    for (size_t i = 0; i < x_bout_interp.size(); ++i) {
        double val = std::sqrt(
            x_bout_interp[i] * x_bout_interp[i] +
            y_bout_interp[i] * y_bout_interp[i] +
            z_bout_interp[i] * z_bout_interp[i]
        );
        vm_bout_interp.push_back(val);
    }

    // Check if mean(vm_bout_interp) > 5 to standardize units
    if (!vm_bout_interp.empty()) {
        double sum_vm = std::accumulate(vm_bout_interp.begin(), vm_bout_interp.end(), 0.0);
        double mean_vm = sum_vm / vm_bout_interp.size();
        if (mean_vm > 5) {
            for (size_t i = 0; i < x_bout_interp.size(); ++i) {
                x_bout_interp[i] /= 9.80665;
                y_bout_interp[i] /= 9.80665;
                z_bout_interp[i] /= 9.80665;
            }
        }
    }

    // Recalculate vm_bout_interp after unit verification and subtract 1
    vm_bout_interp.clear();
    vm_bout_interp.reserve(x_bout_interp.size());
    for (size_t i = 0; i < x_bout_interp.size(); ++i) {
        double val = std::sqrt(
            x_bout_interp[i] * x_bout_interp[i] +
            y_bout_interp[i] * y_bout_interp[i] +
            z_bout_interp[i] * z_bout_interp[i]
        ) - 1.0;
        vm_bout_interp.push_back(val);
    }

    return {t_bout_interp, vm_bout_interp};
}


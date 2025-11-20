#pragma once
#include <vector>
#include <stdexcept>
#include <cmath>
#include <cstddef>
#include <algorithm>
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif


inline std::vector<double> interp1d(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& x_new
) {
     if (x.size() != y.size()) {
        throw std::runtime_error("interp1d: x and y must have the same length"); }
    std::vector<double> y_new(x_new.size());
    
    for (size_t i = 0; i < x_new.size(); ++i) {
        double xi = x_new[i];
        
        // Find the two points to interpolate between
        auto it = std::lower_bound(x.begin(), x.end(), xi);
        
        if (it == x.begin()) {
            y_new[i] = y[0];
        } else if (it == x.end()) {
            y_new[i] = y.back();
        } else {
            size_t idx = std::distance(x.begin(), it);
            double x1 = x[idx - 1];
            double x2 = x[idx];
            double y1 = y[idx - 1];
            double y2 = y[idx];
            
            // Linear interpolation formula
            y_new[i] = y1 + (y2 - y1) * (xi - x1) / (x2 - x1);
        }
    }
    
    return y_new;
}

inline std::vector<double> tukey_window(size_t N, double alpha) {
    std::vector<double> window(N);
    if (N < 2) return window;

    if (alpha <= 0.0) {
        // rectangular
        return window;
    }

    if (alpha >= 1.0) {
        // pure Hann
        for (size_t n = 0; n < N; ++n) {
            window[n] = 0.5 * (1 - std::cos(2 * M_PI * n / (N - 1)));
        }
        return window;
    }

    double alpha_n = alpha * (N - 1) / 2.0;
    
    for (size_t i = 0; i < N; ++i) {
        if (i < alpha_n) {
            window[i] = 0.5 * (1.0 + std::cos(M_PI * (i / alpha_n - 1.0)));
        } else if (i > (N - 1) * (1.0 - alpha / 2.0)) {
            double val = (i - (N - 1) * (1.0 - alpha / 2.0)) / alpha_n;
            window[i] = 0.5 * (1.0 + std::cos(M_PI * val));
        } else {
            window[i] = 1.0;
        }
    }
    
    return window;
}

inline void find_peaks(
    const std::vector<double>& data,
    std::vector<std::size_t>& locs,
    std::vector<double>& pks
) {
    locs.clear();
    pks.clear();

    if (data.size() < 3) {
        return;
    }

    for (std::size_t i = 1; i + 1 < data.size(); ++i) {
        if (data[i] > data[i - 1] && data[i] > data[i + 1]) {
            locs.push_back(i);
            pks.push_back(data[i]);
        }
    }
}
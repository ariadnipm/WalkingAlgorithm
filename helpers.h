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
/*
inline void find_peaks(
    const std::vector<double>& data,
    std::vector<std::size_t>& locs,
    std::vector<double>& pks
) {
    locs.clear();
    pks.clear();

    const std::size_t N = data.size();
    if (N < 3) return;

    std::size_t i = 1;

    while (i + 1 < N) {
       
        if (data[i] > data[i - 1]) {

            
            if (data[i] > data[i + 1]) {
                locs.push_back(i);
                pks.push_back(data[i]);
                ++i;  
            }
            
            else if (data[i] == data[i + 1]) {
                std::size_t left = i;
                std::size_t right = i + 1;

                
                while (right + 1 < N && data[right] == data[right + 1]) {
                    ++right;
                }

                
                if (right < N - 1 && data[right] > data[right + 1]) {
                    std::size_t peak_idx = (left + right) / 2; 
                    locs.push_back(peak_idx);
                    pks.push_back(data[peak_idx]);
                }

                i = right + 1; 
            }
            
            else {
                ++i;
            }
        } else {
            ++i;
        }
    }
} */
inline void find_peaks(
    const std::vector<double>& x,
    std::vector<std::size_t>& locs,
    std::vector<double>& pks
) {
    locs.clear();
    pks.clear();

    const std::size_t N = x.size();
    if (N < 3) return;  // χρειάζονται τουλάχιστον 3 σημεία για peak

    // Υπολογίζουμε "παράγωγο" τύπου np.diff(x)
    std::vector<double> dx(N - 1);
    for (std::size_t i = 0; i + 1 < N; ++i) {
        dx[i] = x[i + 1] - x[i];
    }

    std::size_t i = 1; // αντίστοιχο με το index στο x (1..N-2)

    while (i + 1 < N) {
        // --- Case 1: κανονικό peak (rising then falling) ---
        if (dx[i - 1] > 0.0 && dx[i] < 0.0) {
            std::size_t peak_idx = i;
            locs.push_back(peak_idx);
            pks.push_back(x[peak_idx]);
            ++i;
        }
        // --- Case 2: plateau peak (rising, then flat, μετά falling) ---
        else if (dx[i - 1] > 0.0 && dx[i] == 0.0) {
            std::size_t left = i;
            std::size_t right = i;

            // επέκταση plateau προς τα δεξιά όσο dx == 0
            while (right + 1 < N && dx[right] == 0.0) {
                ++right;
            }
            // τώρα plateau indices στο x είναι [left .. right]
            // dx[right] είναι η πρώτη διαφορά μετά το plateau

            if (right + 1 < N && dx[right] < 0.0) {
                // πρόκειται για plateau peak, πάρε το κέντρο (όπως SciPy)
                std::size_t plateau_left = left;
                std::size_t plateau_right = right;
                std::size_t peak_idx = (plateau_left + plateau_right) / 2;

                locs.push_back(peak_idx);
                pks.push_back(x[peak_idx]);
            }

            // προχώρα μετά το plateau
            i = right + 1;
        }
        else {
            ++i;
        }
    }
}


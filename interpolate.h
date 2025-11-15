#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>

inline std::vector<double> interp1d(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& xnew
) {
    if (x.size() != y.size() || x.size() < 2) {
    
        return {};
    }

    std::vector<double> ynew;
    ynew.reserve(xnew.size());

    size_t j = 0;  

    for (double xn : xnew) {

     
        if (xn <= x.front()) {
            ynew.push_back(y.front());
            continue;
        }

    
        if (xn >= x.back()) {
            ynew.push_back(y.back());
            continue;
        }

        
        while (j + 1 < x.size() && x[j + 1] < xn) {
            ++j;
        }

        double x0 = x[j];
        double x1 = x[j + 1];

        double y0 = y[j];
        double y1 = y[j + 1];

        // linear interpolation
        double t = (xn - x0) / (x1 - x0);
        double yn = y0 + t * (y1 - y0);

        ynew.push_back(yn);
    }

    return ynew;
}
#pragma once
#include <vector>
#include <cmath>

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






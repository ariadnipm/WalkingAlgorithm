#include <iostream>
#include <vector>
#include "cadencefun.h"   
int main() {
    // original x and y values (όπως t_bout και z_bout)
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> y = {0.0, 10.0, 20.0, 30.0};

    // new points where we want interpolation (όπως t_bout_interp)
    std::vector<double> xnew = {0.5, 1.5, 2.5, 3.0, 0.0};

    std::vector<double> ynew = interp1d(x, y, xnew);

    std::cout << "Interpolated values:\n";
    for (size_t i = 0; i < xnew.size(); ++i) {
        std::cout << "x = " << xnew[i] << " → y = " << ynew[i] << "\n";
    }

    return 0;
}
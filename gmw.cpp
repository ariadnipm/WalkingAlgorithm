#include "gmw.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

GeneralizedMorse::GeneralizedMorse(float beta, float gamma)
    : beta_(beta), gamma_(gamma)
{
    if (beta_ <= 0.0f || gamma_ <= 0.0f) {
        std::cerr << "GeneralizedMorse: beta and gamma must be positive; "
                  << "falling back to beta=90, gamma=3\n";
        beta_  = 90.0f;
        gamma_ = 3.0f;
    }

    // modal / peak radian frequency
    float wc = std::pow(beta_ / gamma_, 1.0f / gamma_);
    four_wavelen = wc;          // για Scales 

    imag_frequency = false;     // mother(w) πραγματική στο freq-domain
    doublesided    = false;     // analytic: μόνο θετικά w

    mother = nullptr;
    width  = 0;
}

GeneralizedMorse::~GeneralizedMorse()
{
    if (mother) {
        std::free(mother);
        mother = nullptr;
    }
}

void GeneralizedMorse::generate(int size)
{
    width = size;

    if (mother) {
        std::free(mother);
        mother = nullptr;
    }

    mother = static_cast<float*>(std::malloc(sizeof(float) * size));
    if (!mother) {
        std::cerr << "GeneralizedMorse::generate: malloc failed\n";
        return;
    }

    const float gamma = gamma_;
    const float beta  = beta_;

    const float wc    = std::pow(beta / gamma, 1.0f / gamma);
    const float wcl   = std::log(wc);
    const float wc_gamma = std::pow(wc, gamma);

    // FFT radian step
    const float dw = 2.0f * M_PI / static_cast<float>(size);
    const int   half = size / 2;

    for (int k = 0; k < size; ++k) {
       
        float w;
        if (k <= half) {
            w = k * dw;                      // 0 .. π
        } else {
            w = (k - size) * dw;             // αρνητικές
        }

        if (w <= 0.0f) {
            mother[k] = 0.0f;                // analytic: τίποτα σε αρνητικές & w=0
            continue;
        }

        const float w_gamma = std::pow(w, gamma);

        
        const float exponent = -beta * wcl + wc_gamma +
                               beta * std::log(w) - w_gamma;

        float val = 0.0f;
        if (exponent > -80.0f && exponent < 80.0f) {
            val = 2.0f * std::exp(exponent);
            if (!std::isfinite(val))
                val = 0.0f;
        } else {
            val = 0.0f;
        }

        mother[k] = val;
    }
}

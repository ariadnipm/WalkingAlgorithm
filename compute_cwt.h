#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include "gmw.h"
#include "fCWT/src/fcwt/fcwt.h"
#include "helpers.h"
#include <iomanip>

void print_coefs(const std::vector<std::vector<double>>& C,
                 int max_f = 5, int max_t = 20)
{
    int F = C.size();
    if (F == 0) return;

    int T = C[0].size();

    max_f = std::min(max_f, F);
    max_t = std::min(max_t, T);

    std::cout << "\n=== coefs_interp (partial print) ===\n";
    for (int fi = 0; fi < max_f; ++fi) {
        std::cout << "freq row " << fi << ": ";
        for (int tt = 0; tt < max_t; ++tt) {
            std::cout << std::fixed << std::setprecision(10)
                      << C[fi][tt] << " ";
        }
        std::cout << "...\n";
    }
    std::cout << std::endl;
}

struct CWTResult {
    std::vector<double> freqs_interp;
    std::vector<std::vector<double>> coefs_interp;
};

CWTResult compute_interpolate_cwt(const std::vector<double>& bout,
                                  int fs /*= 10*/)
{
    CWTResult result;

    // αν το bout είναι πολύ μικρό, γύρνα κενό
    if (fs <= 0 || bout.size() < static_cast<std::size_t>(2 * fs)) {
        return result;
    }

    const std::size_t N = bout.size();

    // 1) Tukey window (alpha=0.02 όπως στην Python)
    std::vector<double> window = tukey_window(N, 0.02);
    std::vector<double> windowed_bout(N);
    for (std::size_t i = 0; i < N; ++i) {
        windowed_bout[i] = bout[i] * window[i];
    }

    // 2) zero padding 5*fs αριστερά-δεξιά
    const int padding = 5 * fs;
    std::vector<double> padded_signal;
    padded_signal.reserve(N + 2 * static_cast<std::size_t>(padding));
    padded_signal.insert(padded_signal.end(), padding, 0.0);
    padded_signal.insert(padded_signal.end(),
                         windowed_bout.begin(), windowed_bout.end());
    padded_signal.insert(padded_signal.end(), padding, 0.0);

    // όπως Python: x[:-1]
    if (!padded_signal.empty()) {
        padded_signal.pop_back();
    }
    const int padded_N = static_cast<int>(padded_signal.size());

    // 3) CWT με Generalized Morse (β=90, γ=3)
    const int   n_freqs = 193;
    const float f0      = 0.32f;
    const float f1      = 5.0f;

    GeneralizedMorse gmw(90.0f, 3.0f);
    Scales scales(&gmw, FCWT_LOGSCALES, static_cast<float>(fs),
                  f0, f1, n_freqs);

    // παίρνουμε τις συχνότητες (για παρεμβολή αργότερα)
    std::vector<float> freqs_f(n_freqs);
    scales.getFrequencies(freqs_f.data(), n_freqs);
    std::vector<double> freqs(freqs_f.begin(), freqs_f.end());

   
    std::vector<float> scales_f(n_freqs);
    scales.getScales(scales_f.data(), n_freqs);

    // input -> float
    std::vector<float> sig_f(padded_N);
    for (int i = 0; i < padded_N; ++i) {
        sig_f[i] = static_cast<float>(padded_signal[i]);
    }

    // output: n_freqs x padded_N (scale-major)
    std::vector<std::complex<float>> cwt_out(
        static_cast<std::size_t>(n_freqs) *
        static_cast<std::size_t>(padded_N));

    // use_norm=true: global FFT normalization (ίδιο για όλα τα scales)
    FCWT fcwt_obj(&gmw, /*threads*/1, /*use_opt*/false, /*use_norm*/true);
    fcwt_obj.cwt(sig_f.data(), padded_N, cwt_out.data(), &scales);

    // 4) |CWT|^2 + append last column
    const int Ntime     = padded_N;
    const int Ntime_ext = Ntime + 1;

    std::vector<std::vector<double>> mag2(
        n_freqs, std::vector<double>(Ntime_ext, 0.0));

    for (int f = 0; f < n_freqs; ++f) {
        for (int t = 0; t < Ntime; ++t) {
            const auto& z = cwt_out[static_cast<std::size_t>(f) * Ntime + t];
            const double re = static_cast<double>(z.real());
            const double im = static_cast<double>(z.imag());
            mag2[f][t] = re * re + im * im;  // |W|^2
        }
        mag2[f][Ntime_ext - 1] = mag2[f][Ntime - 1];
    }

   
   
    for (int f = 0; f < n_freqs; ++f) {
        double s = static_cast<double>(scales_f[f]);  // scale από fcwt
        if (s <= 0.0) {
            continue;
        }
        double factor = 1.0 / s;  // αντιστοιχεί σε W / sqrt(s) στο amplitude

        for (int t = 0; t < Ntime_ext; ++t) {
            mag2[f][t] *= factor;
        }
    } 

    // 4c) ensure freqs ascending (αν είναι ανάποδα)
    if (freqs.size() > 1 && freqs[1] < freqs[0]) {
        std::reverse(freqs.begin(), freqs.end());
        std::reverse(mag2.begin(), mag2.end());
        // scales_f δεν το χρειαζόμαστε πια, οπότε δεν πειράζει αν μείνει ως έχει
    }

    // 5) freq grid 0.5–4.45 Hz / 0.05 (όπως Python)
    std::vector<double> freqs_interp;
    const double f_start = 0.5;
    const double f_end   = 4.5;
    const double step    = 0.05;

    const int nfi = static_cast<int>(std::round((f_end - f_start) / step));
    freqs_interp.reserve(nfi);
    for (int i = 0; i < nfi; ++i) {
        freqs_interp.push_back(f_start + step * i);
    }

    // 6) interpolation: coefs_interp[freq_idx][time_idx]
    std::vector<std::vector<double>> coefs_interp(
        nfi, std::vector<double>(Ntime_ext, 0.0));

    for (int t = 0; t < Ntime_ext; ++t) {
        std::vector<double> col(n_freqs);
        for (int f = 0; f < n_freqs; ++f) {
            col[f] = mag2[f][t];
        }

        auto col_interp = interp1d(freqs, col, freqs_interp);
        for (int fi = 0; fi < nfi; ++fi) {
            coefs_interp[fi][t] = col_interp[fi];
        }
    }

    // 7) Trim από COI (5*fs αριστερά–δεξιά)
    const int trim = 5 * fs;
    if (Ntime_ext <= 2 * trim) {
        result.freqs_interp = std::move(freqs_interp);
        result.coefs_interp = std::move(coefs_interp);
        return result;
    }

    const int t_start = trim;
    const int t_end   = Ntime_ext - trim;
    const int Ntrim   = t_end - t_start;

    std::vector<std::vector<double>> coefs_trim(
        nfi, std::vector<double>(Ntrim));

    for (int fi = 0; fi < nfi; ++fi) {
        for (int tt = 0; tt < Ntrim; ++tt) {
            coefs_trim[fi][tt] = coefs_interp[fi][t_start + tt];
        }
    }

    result.freqs_interp = std::move(freqs_interp);
    result.coefs_interp = std::move(coefs_trim);

    std::cout << std::scientific << std::setprecision(6);
   
    return result;
}

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

#include "fCWT/src/fcwt/fcwt.h"
#include "helpers.h"   // interp1d, tukey_window κλπ.

struct CWTResult {
    std::vector<double> freqs_interp;                  // συχνότητες (Hz)
    std::vector<std::vector<double>> coefs_interp;     // [freq][time]
};

CWTResult compute_interpolate_cwt(const std::vector<double>& bout,
                                  int fs = 10)
{
    CWTResult result;

    // αν το bout είναι πολύ μικρό, γύρνα κενά
    if (bout.size() < static_cast<size_t>(2 * fs)) {
        return result;
    }

    const size_t N = bout.size();

    // 1) Tukey window (όπως στην Python: alpha=0.02)
    std::vector<double> window = tukey_window(N, 0.02);
    std::vector<double> windowed_bout(N);
    for (size_t i = 0; i < N; ++i) {
        windowed_bout[i] = bout[i] * window[i];
    }

    // 2) Padding 5*fs zeros σε κάθε πλευρά
    const int padding = 5 * fs;
    std::vector<double> padded_signal;
    padded_signal.reserve(N + 2 * padding);
    padded_signal.insert(padded_signal.end(), padding, 0.0);
    padded_signal.insert(padded_signal.end(),
                         windowed_bout.begin(), windowed_bout.end());
    padded_signal.insert(padded_signal.end(), padding, 0.0);

    // Python: ssq_cwt(tapered_bout[:-1], ...)
    if (!padded_signal.empty()) {
        padded_signal.pop_back();   // αφαιρούμε το τελευταίο sample
    }
    const int padded_N = static_cast<int>(padded_signal.size());

    // 3) CWT με Morlet (fCWT)
    const int   n_freqs = 80;   // ίδιο πλήθος με out[2] στην Python
    const float f0      = 0.5f;
    const float f1      = 4.5f;

    // Morlet wavelet
    Morlet morlet(10.0f);  // fb ≈ 2

    // Scales: LINFREQS για γραμμικό grid 0.5–4.5 Hz
    Scales scales(&morlet, FCWT_LINFREQS, fs, f0, f1, n_freqs);

    // διαβάζουμε τις συχνότητες από τα scales
    std::vector<float> freqs_f(n_freqs);
    scales.getFrequencies(freqs_f.data(), n_freqs);
    std::vector<double> freqs(freqs_f.begin(), freqs_f.end());

    // input σε float για fCWT
    std::vector<float> sig_f(padded_N);
    for (int i = 0; i < padded_N; ++i) {
        sig_f[i] = static_cast<float>(padded_signal[i]);
    }

    // output: n_freqs x padded_N (scale-major)
    std::vector<std::complex<float>> cwt_out(n_freqs * padded_N);

    // FCWT(wavelet*, threads, use_opt_schemes, use_normalization)
    FCWT fcwt_obj(&morlet, /*threads*/1, /*use_opt*/false, /*use_norm*/true);
    fcwt_obj.cwt(sig_f.data(), padded_N, cwt_out.data(), &scales);

    // 4) |CWT|^2 + "append last column" όπως: coefs = np.append(coefs, coefs[:, -1:], 1)
    const int Ntime     = padded_N;
    const int Ntime_ext = Ntime + 1;  // έξτρα στήλη

    std::vector<std::vector<double>> mag2(
        n_freqs, std::vector<double>(Ntime_ext, 0.0)
    );

    for (int f = 0; f < n_freqs; ++f) {
        for (int t = 0; t < Ntime; ++t) {
            const auto& z = cwt_out[f * Ntime + t];
            double re = static_cast<double>(z.real());
            double im = static_cast<double>(z.imag());
            mag2[f][t] = re * re + im * im;  // |z|^2
        }
        // τελευταία στήλη = τελευταία τιμή
        mag2[f][Ntime_ext - 1] = mag2[f][Ntime - 1];
    }

    // 🔁 FIX: αν οι συχνότητες είναι φθίνουσες, γύρνα freqs & mag2
    if (freqs.size() > 1 && freqs[1] < freqs[0]) {
        std::reverse(freqs.begin(), freqs.end());   // 0.5, 0.55, ..., 4.5
        std::reverse(mag2.begin(), mag2.end());     // αντιστοίχιση rows -> freqs
    }

    // 5) freqs_interp = np.arange(0.5, 4.5, 0.05)  (0.5 … 4.45)
    std::vector<double> freqs_interp;
    const double f_start = 0.5;
    const double f_end   = 4.5;
    const double step    = 0.05;

    int nfi = static_cast<int>(std::round((f_end - f_start) / step)); // 80
    freqs_interp.reserve(nfi);
    for (int i = 0; i < nfi; ++i) {
        freqs_interp.push_back(f_start + step * i); // 0.5, 0.55, ..., 4.45
    }

    // coefs_interp[freq_idx][time_idx]
    std::vector<std::vector<double>> coefs_interp(
        nfi, std::vector<double>(Ntime_ext, 0.0)
    );

    // για κάθε χρόνο, interpolάρουμε τη στήλη mag2[:, t] πάνω σε freqs_interp
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

    // 6) Trim από COI: coefs_interp = coefs_interp[:, 5*fs:-5*fs]
    const int trim = 5 * fs;
    if (Ntime_ext <= 2 * trim) {
        // πολύ μικρό → μην κάνεις trim
        result.freqs_interp = std::move(freqs_interp);
        result.coefs_interp = std::move(coefs_interp);
        return result;
    }

    const int t_start = trim;
    const int t_end   = Ntime_ext - trim;  // exclusive
    const int Ntrim   = t_end - t_start;

    std::vector<std::vector<double>> coefs_trim(
        nfi, std::vector<double>(Ntrim)
    );

    for (int fi = 0; fi < nfi; ++fi) {
        for (int tt = 0; tt < Ntrim; ++tt) {
            coefs_trim[fi][tt] = coefs_interp[fi][t_start + tt];
        }
    }

    result.freqs_interp = std::move(freqs_interp);
    result.coefs_interp = std::move(coefs_trim);
    return result;
}

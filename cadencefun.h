#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <functional>
#include <utility>
#include <stdexcept>
#include <optional>
#include "compute_cwt.h"
#include "helpers.h"


inline std::vector<double> adjust_bout(const std::vector<double>& inarray, int fs) {
    std::vector<double> result = inarray;
    size_t remainder = result.size() % fs;
    
    // If data is available for 70% of the last second, pad
    if (remainder >= static_cast<size_t>(0.7 * fs)) {
        size_t padding = fs - remainder;
        double last_val = result.back();
        for (size_t i = 0; i < padding; ++i) {
            result.push_back(last_val);
        }
    } else {
        // Otherwise, trim to full second
        size_t trim_size = (result.size() / fs) * fs;
        result.resize(trim_size);
    }
    
    return result;
}



std::pair<std::vector<double>, std::vector<double>>  preprocess_bout(
    const std::vector<double>& t_bout,
    const std::vector<double>& x_bout,
    const std::vector<double>& y_bout,
    const std::vector<double>& z_bout,
    int fs
) {
    // Check for minimum data
    if (t_bout.size() < 2 || x_bout.size() < 2 || 
        y_bout.size() < 2 || z_bout.size() < 2) {
        return {{}, {}};
    }
    
    // Create interpolated time vector
    std::vector<double> t_bout_interp;
    double t_start = t_bout[0];
    double t_end = t_bout.back();
    double dt = 1.0 / fs;
    
    for (double t = 0; t < (t_end - t_start); t += dt) {
        t_bout_interp.push_back(t + t_start);
    }
    
    if (t_bout_interp.empty()) {
        return {{}, {}};
    }
    
    // Interpolate x, y, z to new time vector
    auto x_bout_interp = interp1d(t_bout, x_bout, t_bout_interp);
    auto y_bout_interp = interp1d(t_bout, y_bout, t_bout_interp);
    auto z_bout_interp = interp1d(t_bout, z_bout, t_bout_interp);
    
    // Adjust bouts to full seconds
    x_bout_interp = adjust_bout(x_bout_interp, fs);
    y_bout_interp = adjust_bout(y_bout_interp, fs);
    z_bout_interp = adjust_bout(z_bout_interp, fs);
    
    // Calculate number of full seconds
    size_t num_seconds = x_bout_interp.size() / fs;
    size_t trim_size = num_seconds * fs;
    
    // Trim time vector
    t_bout_interp.resize(trim_size);
    
    // Decimate time vector to one sample per second
    std::vector<double> t_decimated;
    for (size_t i = 0; i < t_bout_interp.size(); i += fs) {
        t_decimated.push_back(t_bout_interp[i]);
    }
    t_bout_interp = t_decimated;
    
    // Calculate vector magnitude
    std::vector<double> vm_bout_interp(x_bout_interp.size());
    for (size_t i = 0; i < x_bout_interp.size(); ++i) {
        vm_bout_interp[i] = std::sqrt(
            x_bout_interp[i] * x_bout_interp[i] +
            y_bout_interp[i] * y_bout_interp[i] +
            z_bout_interp[i] * z_bout_interp[i]
        );
    }
    
    // Standardize to gravity units (g) if recorded in m/s²
    if (!vm_bout_interp.empty()) {
        double mean_vm = std::accumulate(vm_bout_interp.begin(), 
                                        vm_bout_interp.end(), 0.0) / 
                        vm_bout_interp.size();
        
        if (mean_vm > 5.0) {
            // Convert from m/s² to g
            for (size_t i = 0; i < x_bout_interp.size(); ++i) {
                x_bout_interp[i] /= 9.80665;
                y_bout_interp[i] /= 9.80665;
                z_bout_interp[i] /= 9.80665;
            }
        }
    }
    
    // Recalculate vector magnitude after unit verification and subtract 1g
    for (size_t i = 0; i < x_bout_interp.size(); ++i) {
        vm_bout_interp[i] = std::sqrt(
            x_bout_interp[i] * x_bout_interp[i] +
            y_bout_interp[i] * y_bout_interp[i] +
            z_bout_interp[i] * z_bout_interp[i]
        ) - 1.0;
    }
    
    return {t_bout_interp, vm_bout_interp};
}

std::vector<double> get_pp(const std::vector<double>& vm_bout, int fs)
{
    size_t total = vm_bout.size();
    size_t num_seconds = total / fs;

    std::vector<double> pp(num_seconds);

    for (size_t col = 0; col < num_seconds; ++col) {

        double min_val = vm_bout[col * fs];
        double max_val = vm_bout[col * fs];

        for (int row = 0; row < fs; ++row) {
            size_t idx = col * fs + row;
            double v = vm_bout[idx];

            if (v < min_val) min_val = v;
            if (v > max_val) max_val = v;
        }

        pp[col] = max_val - min_val;
    }

    return pp;
}

std::vector<std::vector<double>> identify_peaks_in_cwt(
    const std::vector<double>& freqs_interp,                       // shape: [num_rows]
    const std::vector<std::vector<double>>& coefs_interp,          // shape: [num_rows][num_cols]
    int fs = 10,
    std::pair<double,double> step_freq = {1.4, 2.3},
    double alpha = 0.6,
    double beta  = 2.5
) {
    std::vector<std::vector<double>> dp;  // θα επιστρέψουμε [num_rows][num_cols2]

    // Έλεγχοι ασφαλείας
    if (coefs_interp.empty() || freqs_interp.empty()) {
        return dp;
    }
    const std::size_t num_rows = coefs_interp.size();
    const std::size_t num_cols = coefs_interp[0].size();
    if (num_rows != freqs_interp.size() || num_cols == 0 || fs <= 0) {
        return dp;
    }

    // βεβαιώσου ότι όλες οι σειρές έχουν ίδιο μέγεθος
    for (const auto& row : coefs_interp) {
        if (row.size() != num_cols) {
            return dp; // μη-ορθογώνιος πίνακας => error
        }
    }

    const std::size_t num_cols2 = num_cols / static_cast<std::size_t>(fs);
    if (num_cols2 == 0) {
        return dp;
    }

    dp.assign(num_rows, std::vector<double>(num_cols2, 0.0));

    // εύρεση loc_min, loc_max όπως: np.argmin(abs(freqs_interp - step_freq[0]))
    auto argmin_abs = [](const std::vector<double>& v, double target) {
        std::size_t idx_min = 0;
        double best = std::fabs(v[0] - target);
        for (std::size_t i = 1; i < v.size(); ++i) {
            double diff = std::fabs(v[i] - target);
            if (diff < best) {
                best = diff;
                idx_min = i;
            }
        }
        return idx_min;
    };

    const std::size_t loc_min = argmin_abs(freqs_interp, step_freq.first);
    const std::size_t loc_max = argmin_abs(freqs_interp, step_freq.second);

    // Loop σε κάθε 1-second μη-επικαλυπτόμενο παράθυρο στον χρόνο
    for (std::size_t i = 0; i < num_cols2; ++i) {
        const std::size_t x_start = i * static_cast<std::size_t>(fs);
        const std::size_t x_end   = x_start + static_cast<std::size_t>(fs);

        if (x_end > num_cols) {
            break; // ασφαλείας
        }

        // window = sum(coefs_interp[:, x_start:x_end], axis=1)
        std::vector<double> window(num_rows, 0.0);
        for (std::size_t r = 0; r < num_rows; ++r) {
            double sum = 0.0;
            for (std::size_t t = x_start; t < x_end; ++t) {
                sum += coefs_interp[r][t];
            }
            window[r] = sum;
        }

        // find_peaks(window) -> locs, pks
        std::vector<std::size_t> locs;
        std::vector<double> pks;
        find_peaks(window, locs, pks);

        if (locs.empty()) {
            // δεν υπάρχουν peaks σ' αυτό το παράθυρο
            continue;
        }

        // ταξινόμηση peaks κατά φθίνουσα τιμή, όπως:
        // ind = np.argsort(-pks)
        std::vector<std::size_t> order(locs.size());
        for (std::size_t k = 0; k < order.size(); ++k) {
            order[k] = k;
        }

        std::sort(order.begin(), order.end(),
                  [&](std::size_t a, std::size_t b) {
                      return pks[a] > pks[b]; // descending
                  });

        std::vector<std::size_t> locs_sorted(locs.size());
        std::vector<double>      pks_sorted(pks.size());
        for (std::size_t k = 0; k < order.size(); ++k) {
            locs_sorted[k] = locs[order[k]];
            pks_sorted[k]  = pks[order[k]];
        }
        locs.swap(locs_sorted);
        pks.swap(pks_sorted);

        // βρες πρώτο peak μέσα στο [loc_min, loc_max]
        std::optional<std::size_t> index_in_range_opt;
        for (std::size_t j = 0; j < locs.size(); ++j) {
            std::size_t loc_j = locs[j];
            if (loc_j >= loc_min && loc_j <= loc_max) {
                index_in_range_opt = j;
                break;
            }
        }

        std::vector<double> peak_vec(num_rows, 0.0);

        if (index_in_range_opt.has_value()) {
            std::size_t index_in_range = *index_in_range_opt;

            // peaks_below = pks[locs < loc_min]
            // peaks_above = pks[locs > loc_max]
            double max_peak_below = 0.0;
            double max_peak_above = 0.0;

            for (std::size_t j = 0; j < locs.size(); ++j) {
                std::size_t loc_j = locs[j];
                double pk = pks[j];
                if (loc_j < loc_min) {
                    if (pk > max_peak_below) max_peak_below = pk;
                } else if (loc_j > loc_max) {
                    if (pk > max_peak_above) max_peak_above = pk;
                }
            }

            double center_peak = pks[index_in_range];

            if (center_peak > 0.0) {
                double ratio_b = max_peak_above / center_peak;
                double ratio_a = max_peak_below / center_peak;

                // if max_b / center < beta OR max_a / center < alpha:
                if (ratio_b < beta || ratio_a < alpha) {
                    std::size_t freq_idx = locs[index_in_range];
                    if (freq_idx < num_rows) {
                        peak_vec[freq_idx] = 1.0;
                    }
                }
            }
        }

        // dp[:, i] = peak_vec
        for (std::size_t r = 0; r < num_rows; ++r) {
            dp[r][i] = peak_vec[r];
        }
    }

    return dp;
}







std::vector<std::vector<double>> find_continuous_dominant_peaks(
    const std::vector<std::vector<double>>& valid_peaks,
    int min_t,
    int delta
) {
    std::vector<std::vector<double>> empty;  // για early return

    if (valid_peaks.empty() || min_t <= 0) {
        return empty;
    }

    const std::size_t num_rows = valid_peaks.size();
    const std::size_t num_cols = valid_peaks[0].size();

    if (num_cols == 0) {
        return empty;
    }

    // έλεγχος ότι όλες οι σειρές έχουν ίδιο μήκος
    for (const auto& row : valid_peaks) {
        if (row.size() != num_cols) {
            return empty;  // μη ορθογώνιος πίνακας
        }
    }

    // extended_peaks: [num_rows][num_cols+1]
    std::vector<std::vector<double>> extended_peaks(
        num_rows, std::vector<double>(num_cols + 1, 0.0)
    );

    for (std::size_t r = 0; r < num_rows; ++r) {
        for (std::size_t c = 0; c < num_cols; ++c) {
            extended_peaks[r][c] = valid_peaks[r][c];
        }
    }

    // cont_peaks: ίδιο shape με extended
    std::vector<std::vector<double>> cont_peaks(
        num_rows, std::vector<double>(num_cols + 1, 0.0)
    );

    const int num_cols_ext = static_cast<int>(num_cols) + 1;

    // windows = [0..min_t-1] + [min_t-2..0]
    std::vector<int> windows;
    windows.reserve(2 * min_t - 1);
    for (int w = 0; w < min_t; ++w) {
        windows.push_back(w);
    }
    for (int w = min_t - 2; w >= 0; --w) {
        windows.push_back(w);
    }

    // κύριος βρόχος πάνω στα slice indices
    // range(num_cols + 1 - min_t)
    for (int slice_ind = 0; slice_ind <= num_cols_ext - min_t - 1; ++slice_ind) {

        // slice_mat: [num_rows][min_t]
        std::vector<std::vector<double>> slice_mat(
            num_rows, std::vector<double>(min_t, 0.0)
        );

        for (std::size_t r = 0; r < num_rows; ++r) {
            for (int t = 0; t < min_t; ++t) {
                slice_mat[r][t] = extended_peaks[r][slice_ind + t];
            }
        }

        bool stop = true;

        // for win_ind in windows:
        for (int w_idx = 0; w_idx < static_cast<int>(windows.size()); ++w_idx) {
            int win_ind = windows[w_idx];

            stop = true;  // όπως στην Python, ξαναγίνεται True στην αρχή κάθε iteration

            // pr = where(slice_mat[:, win_ind] != 0)
            std::vector<int> pr;
            for (std::size_t r = 0; r < num_rows; ++r) {
                if (slice_mat[r][win_ind] != 0.0) {
                    pr.push_back(static_cast<int>(r));
                }
            }

            // αν δεν υπάρχουν peaks σε αυτή τη στήλη
            if (pr.empty()) {
                // stop μένει True -> break όπως στην Python
                break;
            }

            // επεξεργασία για κάθε peak p
            for (int p_idx = 0; p_idx < static_cast<int>(pr.size()); ++p_idx) {
                int p = pr[p_idx];

                int start_idx = std::max(0, p - delta);
                int end_idx   = std::min(p + delta + 1, static_cast<int>(num_rows));

                int index_len = end_idx - start_idx;
                if (index_len <= 0) {
                    continue;
                }

                double center_val = slice_mat[p][win_ind];

                // peaks1, peaks2 σαν τα NumPy vectors
                std::vector<double> peaks1(index_len, center_val);
                std::vector<double> peaks2(index_len, center_val);

                if (win_ind == 0) {
                    // peaks1 += slice_mat[index, win_ind + 1]
                    for (int k = 0; k < index_len; ++k) {
                        int rr = start_idx + k;
                        peaks1[k] += slice_mat[rr][win_ind + 1];
                    }
                } else if (win_ind == min_t - 1) {
                    // peaks1 += slice_mat[index, win_ind - 1]
                    for (int k = 0; k < index_len; ++k) {
                        int rr = start_idx + k;
                        peaks1[k] += slice_mat[rr][win_ind - 1];
                    }
                } else {
                    // middle case
                    // peaks1 += slice_mat[index, win_ind - 1]
                    // peaks2 += slice_mat[index, win_ind + 1]
                    for (int k = 0; k < index_len; ++k) {
                        int rr = start_idx + k;
                        peaks1[k] += slice_mat[rr][win_ind - 1];
                        peaks2[k] += slice_mat[rr][win_ind + 1];
                    }
                }

                // any(peaks1 > 1) / any(peaks2 > 1)
                auto any_greater_than_one = [](const std::vector<double>& v) {
                    for (double val : v) {
                        if (val > 1.0) return true;
                    }
                    return false;
                };

                bool any1 = any_greater_than_one(peaks1);
                bool any2 = any_greater_than_one(peaks2);

                if (win_ind == 0 || win_ind == min_t - 1) {
                    if (any1) {
                        stop = false;
                    } else {
                        slice_mat[p][win_ind] = 0.0;
                    }
                } else {
                    if (any1 && any2) {
                        stop = false;
                    } else {
                        slice_mat[p][win_ind] = 0.0;
                    }
                }
            } // end for p in pr

            if (stop) {
                // όπως στο Python: if stop: break
                break;
            }

        } // end for win_ind in windows

        if (!stop) {
            // cont_peaks[:, slice_ind:slice_ind+min_t] = slice_mat
            for (std::size_t r = 0; r < num_rows; ++r) {
                for (int t = 0; t < min_t; ++t) {
                    cont_peaks[r][slice_ind + t] = slice_mat[r][t];
                }
            }
        }
    }

    // return cont_peaks[:, :-1]
    std::vector<std::vector<double>> cont_trim(
        num_rows, std::vector<double>(num_cols, 0.0)
    );

    for (std::size_t r = 0; r < num_rows; ++r) {
        for (std::size_t c = 0; c < num_cols; ++c) {
            cont_trim[r][c] = cont_peaks[r][c];
        }
    }

    return cont_trim;
}



std::vector<double> find_walking(
    const std::vector<double>& vm_bout,
    int fs        = 10,
    double min_amp = 0.3,
    std::pair<double,double> step_freq = {1.4, 2.3},
    double alpha  = 0.6,
    double beta   = 2.5,
    int min_t     = 3,
    int delta     = 20
) {
    // πόσα "δευτερόλεπτα" έχουμε (υποθέτουμε ότι vm_bout.length % fs == 0)
    const std::size_t n_seconds = vm_bout.size() / fs;

    // --- 1) peak-to-peak ανά δευτερόλεπτο (pp) ---
    std::vector<double> pp = get_pp(vm_bout, fs);   // μέγεθος n_seconds

    // --- 2) μάσκα valid (high-intensity) ---
    std::vector<char> valid(n_seconds, 1);          // 1 = True, 0 = False
    for (std::size_t i = 0; i < n_seconds; ++i) {
        if (pp[i] < min_amp) {
            valid[i] = 0;
        }
    }

    int sum_valid = std::accumulate(valid.begin(), valid.end(), 0);

    // αποτέλεσμα: cadence (Hz) ανά δευτερόλεπτο
    std::vector<double> cad(n_seconds, 0.0);

    // --- 3) αν δεν έχουμε αρκετά valid seconds, επιστρέφουμε μηδενικά ---
    if (sum_valid < min_t) {
        // cad είναι ήδη γεμάτο με 0
        return cad;
    }

    // --- 4) φτιάχνουμε το tapered_bout = vm_bout[np.repeat(valid, fs)] ---
    std::vector<double> tapered_bout;
    tapered_bout.reserve(static_cast<std::size_t>(sum_valid) * fs);

    for (std::size_t sec = 0; sec < n_seconds; ++sec) {
        if (valid[sec]) {
            for (int k = 0; k < fs; ++k) {
                std::size_t idx = sec * fs + k;
                if (idx < vm_bout.size()) {
                    tapered_bout.push_back(vm_bout[idx]);
                }
            }
        }
    }

    // --- 5) CWT + interpolation ---
    CWTResult cwt = compute_interpolate_cwt(tapered_bout, fs);
    const auto& freqs_interp = cwt.freqs_interp;            // μέγεθος n_freqs
    const auto& coefs_interp = cwt.coefs_interp;            // [freq][time]

    if (freqs_interp.empty() || coefs_interp.empty()) {
        return cad; // όλα 0
    }

    // --- 6) Dominant peaks στον CWT ---
    auto dp = identify_peaks_in_cwt(freqs_interp, coefs_interp,
                                    fs, step_freq, alpha, beta);
    // dp: [n_freqs][n_valid_seconds_μετά_το_CWT]
    std::size_t n_freqs = dp.size();
    std::size_t n_dp_cols = dp[0].size();   // # seconds που περιγράφει το dp

    // --- 7) valid_peaks: μοιράζουμε τα dp στις σωστές χρονικές θέσεις ---
    // valid_peaks: [n_freqs][n_seconds]  (ίδιο με Python: (dp.shape[0], len(valid)))
    std::vector<std::vector<double>> valid_peaks(
        n_freqs, std::vector<double>(n_seconds, 0.0)
    );

    // dp αντιστοιχεί ΜΟΝΟ στα "True" seconds, κατά σειρά.
    std::size_t dp_col = 0;
    for (std::size_t sec = 0; sec < n_seconds && dp_col < n_dp_cols; ++sec) {
        if (valid[sec]) {
            for (std::size_t f = 0; f < n_freqs; ++f) {
                valid_peaks[f][sec] = dp[f][dp_col];
            }
            ++dp_col;
        }
    }

    // --- 8) continuous peaks σε χρόνο & συχνότητα ---
    auto cont_peaks = find_continuous_dominant_peaks(valid_peaks,
                                                     min_t, delta);
    // cont_peaks: ίδιο shape με valid_peaks

    // --- 9) Περίληψη: για κάθε δευτερόλεπτο, αν υπάρχει κάποιο continuous peak
    //                  παίρνουμε την πρώτη συχνότητα και τη βάζουμε στο cad
    for (std::size_t sec = 0; sec < n_seconds; ++sec) {
        for (std::size_t f = 0; f < n_freqs; ++f) {
            if (cont_peaks[f][sec] > 0.0) {
                cad[sec] = freqs_interp[f];
                break;  // μόνο την πρώτη (χαμηλότερη) συχνότητα κρατάμε
            }
        }
    }

    return cad;
}
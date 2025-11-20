#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "compute_cwt.h"
#include "cadencefun.h"

  /*
struct TestBout {
    std::vector<double> t;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
};

inline TestBout generate_test_bout(int duration_sec = 10, int fs = 10) {
    

    int N = duration_sec * fs;
    std::vector<double> t(N), x(N), y(N), z(N);

    double f_step = 3.0; // Hz

    for (int i = 0; i < N; ++i) {
        double n = static_cast<double>(i);
        double t_ideal = n / static_cast<double>(fs);

        // timestamps με jitter (ίδιος τύπος με Python)
        t[i] = 1000.0 + t_ideal + 0.003 * std::sin(0.37 * n);

        // σήματα σε g, ίδιοι τύποι με Python:
        x[i] =
            0.3 * std::sin(2.0 * M_PI * f_step * t_ideal) +
            0.05 * std::sin(2.0 * M_PI * 0.5 * t_ideal);

        y[i] =
            0.2 * std::sin(2.0 * M_PI * f_step * t_ideal + 0.7) +
            0.03 * std::sin(2.0 * M_PI * 0.8 * t_ideal);

        z[i] =
            1.0 +
            0.6 * std::sin(2.0 * M_PI * f_step * t_ideal + 0.3);
    }

    return {t, x, y, z};
}
void print_vec(const std::string& name, const std::vector<double>& v) {
    std::cout << name << " (size=" << v.size() << "): ";
    for (double x : v) std::cout << x << " ";
    std::cout << "\n";
}

int main() {
    const int fs = 10;

    // 1) Φτιάχνουμε το τεστ σήμα
    TestBout bout = generate_test_bout(10, fs);

    // 2) Καλούμε preprocess_bout
    auto prep = preprocess_bout(bout.t, bout.x, bout.y, bout.z, fs);
    std::vector<double>& t_final  = prep.first;
    std::vector<double>& vm_final = prep.second;

    // 3) Καλούμε get_pp
    std::vector<double> pp = get_pp(vm_final, fs);
   
   auto result = compute_interpolate_cwt(vm_final, 10);
   auto& freqs_interp = result.freqs_interp;
   auto& coefs_interp = result.coefs_interp;


    // 4) Τυπώνουμε
    print_vec("t_final", t_final);
    print_vec("vm_final", vm_final);
    print_vec("pp", pp);
    print_vec("freqs_interp", freqs_interp);
    print_vec("coefs_interp[0]", coefs_interp[0]); // πρώτη σειράs

std::vector<double> energy_per_freq(coefs_interp.size(), 0.0);
for (size_t fi = 0; fi < coefs_interp.size(); ++fi) {
    for (double v : coefs_interp[fi]) energy_per_freq[fi] += v;
}

size_t argmax = std::distance(
    energy_per_freq.begin(),
    std::max_element(energy_per_freq.begin(), energy_per_freq.end())
);

std::cout << "Peak freq C++: " << freqs_interp[argmax] << " Hz\n";


    return 0;
} */



// Generate synthetic signal: sum of sinusoids

constexpr int FS = 10;
constexpr int DURATION = 60;
constexpr int N = FS * DURATION;

std::vector<double> signal_test_A()
{
    // Test A: δύο συχνότητες 1.8 Hz + 2.6 Hz
    double f1 = 1.8;
    double f2 = 2.6;
    std::vector<double> x(N);
    for (int n = 0; n < N; ++n) {
        double t = static_cast<double>(n) / FS;
        x[n] = std::sin(2.0 * M_PI * f1 * t)
             +   std::sin(2.0 * M_PI * f2 * t);
    }
    return x;
}

std::vector<double> signal_test_B()
{
    // Test B: τρία κομμάτια με διαφορετική συχνότητα
    double f1 = 1.4;
    double f2 = 1.8;
    double f3 = 2.2;

    std::vector<double> x(N);
    int n1 = N / 3;
    int n2 = 2 * N / 3;

    for (int n = 0; n < N; ++n) {
        double t = static_cast<double>(n) / FS;
        if (n < n1) {
            x[n] = std::sin(2.0 * M_PI * f1 * t);
        } else if (n < n2) {
            x[n] = std::sin(2.0 * M_PI * f2 * t);
        } else {
            x[n] = std::sin(2.0 * M_PI * f3 * t);
        }
    }
    return x;
}

std::vector<double> signal_test_C()
{
    // Test C: "ψευδο-τυχαίος" θόρυβος, deterministic (fixed seed-like)
    std::vector<double> x(N);
    // απλός LCG για να μην μπλέκουμε με <random>
    unsigned int seed = 123456789u;
    auto lcg = [&seed]() {
        seed = 1664525u * seed + 1013904223u;
        return seed;
    };

    for (int n = 0; n < N; ++n) {
        // uniform σε (0,1)
        double u = static_cast<double>(lcg()) / static_cast<double>(UINT32_MAX);
        // approx zero-mean "noise"
        x[n] = 0.5 * (u - 0.5);
    }
    return x;
}

void run_test_cpp(const std::vector<double>& x, int fs, const std::string& label)
{
    CWTResult cwt = compute_interpolate_cwt(x, fs);

    const auto& freqs_interp = cwt.freqs_interp;
    const auto& coefs_interp = cwt.coefs_interp;  // [freq][time]

    if (freqs_interp.empty() || coefs_interp.empty()) {
        std::cout << "[C++] " << label << " => EMPTY CWT\n";
        return;
    }

    // ενέργεια ανά συχνότητα
    std::vector<double> energy_per_freq(coefs_interp.size(), 0.0);
    for (size_t fi = 0; fi < coefs_interp.size(); ++fi) {
        for (double v : coefs_interp[fi]) {
            energy_per_freq[fi] += v;
        }
    }

    auto it_max = std::max_element(energy_per_freq.begin(), energy_per_freq.end());
    size_t argmax = std::distance(energy_per_freq.begin(), it_max);
    double peak_freq = freqs_interp[argmax];

    std::cout << "[C++] " << label << " peak freq = " << peak_freq << " Hz\n";

    // προαιρετικά: top-3 peaks
    std::vector<size_t> idx(coefs_interp.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
        return energy_per_freq[a] > energy_per_freq[b];
    });

    std::cout << "   top3 freqs: ";
    for (int k = 0; k < 3 && k < (int)idx.size(); ++k) {
        std::cout << freqs_interp[idx[k]] << " ";
    }
    std::cout << "\n";
}
void print_dp(const std::vector<std::vector<double>>& dp, const std::string& name) {
    std::cout << name << " (rows=" << dp.size() << "):\n";
    for (std::size_t r = 0; r < dp.size(); ++r) {
        std::cout << "row " << r << ": ";
        for (std::size_t c = 0; c < dp[r].size(); ++c) {
            std::cout << dp[r][c] << " ";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

void print_matrix(const std::vector<std::vector<double>>& M, const std::string& name){
    std::cout << name << ":\n";
    for (const auto& row : M){
        for (double x : row) std::cout << x << " ";
        std::cout << "\n";
    }
    std::cout << "\n";
}




int main()
{
   /* auto xA = signal_test_A();
    auto xB = signal_test_B();
    auto xC = signal_test_C();

    run_test_cpp(xA, FS, "Test A (1.8 + 2.6 Hz)");
    run_test_cpp(xB, FS, "Test B (1.4 -> 1.8 -> 2.2 Hz)");
    run_test_cpp(xC, FS, "Test C (noise)");
 
std::vector<double> freqs_interp = {1.0, 1.5, 2.0, 3.0};

std::vector<std::vector<double>> coefs_interp(4, std::vector<double>(20, 0.0));

// 1ο δευτερόλεπτο (0..9)
for (int t = 0; t < 10; ++t) {
    coefs_interp[0][t] = 1.0;  // 1.0 Hz
    coefs_interp[1][t] = 5.0;  // 1.5 Hz (peak)
    coefs_interp[2][t] = 2.0;  // 2.0 Hz
    coefs_interp[3][t] = 1.5;  // 3.0 Hz
}

// 2ο δευτερόλεπτο (10..19)
for (int t = 10; t < 20; ++t) {
    coefs_interp[0][t] = 2.0;  // 1.0 Hz
    coefs_interp[1][t] = 1.0;  // 1.5 Hz
    coefs_interp[2][t] = 6.0;  // 2.0 Hz (peak)
    coefs_interp[3][t] = 5.0;  // 3.0 Hz
}
auto dp = identify_peaks_in_cwt(freqs_interp, coefs_interp);


print_dp(dp, "dp");
*/

   std::vector<std::vector<double>> valid_A(6, std::vector<double>(10, 0.0));
    valid_A[2][3] = 1;
    valid_A[2][4] = 1;
    valid_A[2][5] = 1;

    auto outA = find_continuous_dominant_peaks(valid_A, 3, 1);
    print_matrix(outA, "C++ Test A");

    // === Case B ===
    std::vector<std::vector<double>> valid_B(6, std::vector<double>(10, 0.0));
    valid_B[3][1] = 1;
    valid_B[3][4] = 1;

    auto outB = find_continuous_dominant_peaks(valid_B, 3, 1);
    print_matrix(outB, "C++ Test B");

    // === Case C ===
    std::vector<std::vector<double>> valid_C(6, std::vector<double>(10, 0.0));
    valid_C[1][2] = 1;
    valid_C[1][4] = 1;

    auto outC = find_continuous_dominant_peaks(valid_C, 3, 1);
    print_matrix(outC, "C++ Test C");

    return 0;
}
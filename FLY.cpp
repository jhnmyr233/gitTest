#include <vector>
#include <complex>
#include <cmath>

using std::vector;
using std::complex;

// MATLAB-like FFT / IFFT 接口
vector<complex<double>> fft(const vector<complex<double>>& x, int n = -1) {
    int N = (n > 0 ? n : static_cast<int>(x.size()));
    vector<complex<double>> X(N);
    const double PI = std::acos(-1.0);
    for (int k = 0; k < N; ++k) {
        complex<double> sum(0.0, 0.0);
        for (int m = 0; m < N; ++m) {
            complex<double> value = (m < static_cast<int>(x.size()) ? x[m] : complex<double>(0.0, 0.0));
            double angle = -2.0 * PI * k * m / N;
            sum += value * complex<double>(std::cos(angle), std::sin(angle));
        }
        X[k] = sum;
    }
    return X;
}

vector<complex<double>> fft(const vector<double>& x, int n = -1) {
    vector<complex<double>> xc(x.begin(), x.end());
    return fft(xc, n);
}

vector<complex<double>> ifft(const vector<complex<double>>& X, int n = -1) {
    int N = (n > 0 ? n : static_cast<int>(X.size()));
    vector<complex<double>> x(N);
    const double PI = std::acos(-1.0);
    for (int m = 0; m < N; ++m) {
        complex<double> sum(0.0, 0.0);
        for (int k = 0; k < N; ++k) {
            complex<double> value = (k < static_cast<int>(X.size()) ? X[k] : complex<double>(0.0, 0.0));
            double angle = 2.0 * PI * k * m / N;
            sum += value * complex<double>(std::cos(angle), std::sin(angle));
        }
        x[m] = sum / static_cast<double>(N);
    }
    return x;
}

// MATLAB 风格接口示例：
// vector<complex<double>> Y = fft(x);
// vector<complex<double>> y = ifft(Y);

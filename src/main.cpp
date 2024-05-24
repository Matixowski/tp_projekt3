#define _USE_MATH_DEFINES

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <matplot/matplot.h>
#include <cmath>
#include <vector>
#include <complex>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

std::vector<std::complex<double>> dft(std::vector<double> a, int n) {
    std::complex<double> tmp;
    std::complex<double> z;
    std::complex<double> im(0,1);
    std::vector<std::complex<double>> y;
    for (int i = 0; i < n; i++) {
        z.real(0);
        z.imag(0);
        for (int j = 0; j < n; j++) {
            tmp.real(a[j]);
            tmp.imag(0);
            z += tmp * exp(-2.0 * im * M_PI * (((double)i * j) / n));
        }
        y.push_back(z);
    }
    return y;
}

std::vector<std::complex<double>> idft(std::vector<std::complex<double>> a, int n) {
    std::complex<double> z;
    std::complex<double> im(0, 1);
    std::vector<std::complex<double>> y;
    for (int i = 0; i < n; i++) {
        z.real(0);
        z.imag(0);
        for (int j = 0; j < n; j++) {
            z += a[j] * exp(2.0 * im * M_PI * (((double)i * j) / n));
        }
        z /= (double)n;
        y.push_back(z);
    }
    return y;
}

std::vector<float> deriv(std::vector<float> a, int n) {
    std::vector<float> y;
    for (int i = 0; i < n - 1; i++) {
        y.push_back(a[i+1] - a[i]);
    }
    return y;
}

void spectrum(std::vector<std::complex<double>> a, int n, int samprate, std::string range) {
    std::vector<double> freq;
    std::vector<double> ampl;
    int m = n / 2;
    if (range == "full") {
        m = n;
    }
    for (int i = 0; i < m; i++) {
        freq.push_back((double)samprate * i / n);
        ampl.push_back(sqrt(pow(a[i].real(), 2) + pow(a[i].imag(), 2))/n);
    }
    matplot::plot(freq, ampl);
    matplot::show();
}

void plot(std::vector<double> signal, std::vector<double> t = {}) {
    if (t.size() == 0) {
        matplot::plot(signal);
    }
    else {
        matplot::plot(t, signal);
    }
    matplot::show();
}

double genSquareWave(double t, double frequency, double amplitude, double phase_in_degrees) {
    return amplitude * copysign(1.0, sin(2 * M_PI * (t *frequency - phase_in_degrees / 360)));
}

double genSawtoothWave(double t, double frequency, double amplitude, double phase_in_degrees) {
    return amplitude * 2 / M_PI * atan(tan(M_PI * (t *frequency - (phase_in_degrees + 180) / 360)));
}

double genSineWave(double t, double frequency, double amplitude, double phase_in_degrees) {
    return amplitude * sin(2 * M_PI * (t *frequency - phase_in_degrees / 360));
}

double genSignal(double t, std::string type, double frequency, double amplitude = 1, double phase_in_degrees = 0) {
    if (type == "square") {
        return genSquareWave(t, frequency, amplitude, phase_in_degrees);
    }
    if (type == "sawtooth") {
        return genSawtoothWave(t, frequency, amplitude, phase_in_degrees);
    }
    if (type == "sine") {
        return genSineWave(t, frequency, amplitude, phase_in_degrees);
    }
    if (type == "cosine") {
        phase_in_degrees -= 90;
        return genSineWave(t, frequency, amplitude, phase_in_degrees);
    }
}

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
           dft
           rdft
    )pbdoc";

    m.def("dft", &dft, R"pbdoc(
        Direct Fourier Transform
    )pbdoc");

    m.def("idft", &idft, R"pbdoc(
        Inverse Fourirer Transform
    )pbdoc");

    m.def("spectrum", &spectrum, py::arg(), py::arg(), py::arg(), py::arg("range") = "half", R"pbdoc(
        Display frequency spectrum from DFT
    )pbdoc");

    m.def("deriv", &deriv, R"pbdoc(
        Signal derivative
    )pbdoc");

    std::vector<double> v = {};
    m.def("plot", &plot, py::arg("signal"), py::arg("t") = v, R"pbdoc(
        Plot signal
    )pbdoc");

    m.def("genSignal", py::vectorize(genSignal), py::arg("t"), py::arg("type"), py::arg("frequency"), py::arg("amplitude") = 1, py::arg("phase_in_degrees") = 0, R"pbdoc(
        Generate signal
    )pbdoc");

/*#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif*/
    m.attr("__version__") = "dev";
}

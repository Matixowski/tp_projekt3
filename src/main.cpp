#define _USE_MATH_DEFINES

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <matplot/matplot.h>
#include <cmath>
#include <vector>
#include <complex>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

struct complexnumber {
    double real = 0.0;
    double imaginary = 0.0;
};

int add(int i, int j) {
    return i + j;
}

std::vector<std::complex<double>> dft(py::list a, int n) {
    std::complex<double> tmp;
    std::complex<double> z;
    std::complex<double> im(0,1);
    std::vector<std::complex<double>> y;
    for (int i = 0; i < n; i++) {
        z.real(0);
        z.imag(0);
        for (int j = 0; j < n; j++) {
            tmp.real(a[j].cast<double>());
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

std::vector<int> deriv(py::list a, int n) {
    std::vector<int> y;
    for (int i = 0; i < n - 1; i++) {
        y.push_back(a[i+1].cast<int>() - a[i].cast<int>());
    }
    return y;
}

void spectrum(std::vector<std::complex<double>> a, int n, int samprate) {
    std::vector<double> freq;
    std::vector<double> ampl;
    for (int i = 0; i < n; i++) {
        freq.push_back((double)samprate * i / n);
        ampl.push_back(sqrt(pow(a[i].real(), 2) + pow(a[i].imag(), 2))/n);
    }
    matplot::plot(freq, ampl);
    matplot::show();
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

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.def("dft", &dft, R"pbdoc(
        Direct Fourier Transform
    )pbdoc");

    m.def("idft", &idft, R"pbdoc(
        Inverse Fourirer Transform
    )pbdoc");

    m.def("spectrum", &spectrum, R"pbdoc(
        Display frequency spectrum from DFT
    )pbdoc");

    m.def("deriv", &deriv, R"pbdoc(
        Signal derivative
    )pbdoc");

/*#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif*/
    m.attr("__version__") = "dev";
}

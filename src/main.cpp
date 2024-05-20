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

void plot(std::vector<double> signal) {
    matplot::plot(signal);
    matplot::show();
}

std::vector<double> genSquareWave(double amplitude, int wavelength, int samples, int phase_in_degrees) {
    std::vector<double> wave;
    wave.resize(samples);
    int phase_in_samples = lround(phase_in_degrees / 360.0 * wavelength);
    for (int i = 0; i < wavelength; i++) {
        if ((i + wavelength - phase_in_samples) * 2 / wavelength % 2) {
            wave[i] = -amplitude;
        }
        else {
            wave[i] = amplitude;
        }
    }
    for (int i = wavelength; i < samples; i++) {
        wave[i] = wave[i % wavelength];
    }
    return wave;
}

std::vector<double> genSawtoothWave(double amplitude, int wavelength, int samples, int phase_in_degrees) {
    std::vector<double> wave;
    wave.resize(samples);
    int phase_in_samples = lround(phase_in_degrees / 360.0 * wavelength);
    double step = amplitude * 2 / (wavelength - 1);
    for (int i = 0; i < wavelength; i++) {
        int j = (i + wavelength - phase_in_samples) % wavelength;
        if (j == wavelength - 1) {
            wave[i] = amplitude;
        }
        else {
            wave[i] = -amplitude + j * step;
        }
    }
    for (int i = wavelength; i < samples; i++) {
        wave[i] = wave[i % wavelength];
    }
    return wave;
}

std::vector<double> genSineWave(double amplitude, int wavelength, int samples, int phase_in_degrees) {
    std::vector<double> wave;
    wave.resize(samples);
    for (int i = 0; i < wavelength; i++) {
        wave[i] = amplitude * sin(2 * M_PI * (1.0 * i / wavelength - phase_in_degrees / 360.0));
    }
    for (int i = wavelength; i < samples; i++) {
        wave[i] = wave[i % wavelength];
    }
    return wave;
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

    m.def("spectrum", &spectrum, py::arg(), py::arg(), py::arg(), py::arg("range") = "half", R"pbdoc(
        Display frequency spectrum from DFT
    )pbdoc");

    m.def("deriv", &deriv, R"pbdoc(
        Signal derivative
    )pbdoc");

    m.def("plot", &plot, R"pbdoc(
        Plot signal
    )pbdoc");

    m.def("genSquareWave", &genSquareWave, R"pbdoc(
        Generate square wave
    )pbdoc");

    m.def("genSawtoothWave", &genSawtoothWave, R"pbdoc(
        Generate sawtooth wave
    )pbdoc");

    m.def("genSineWave", &genSineWave, R"pbdoc(
        Generate sine wave
    )pbdoc");

/*#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif*/
    m.attr("__version__") = "dev";
}

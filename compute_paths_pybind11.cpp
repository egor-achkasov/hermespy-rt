#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <complex>
#include <vector>
#include <iostream>

#ifdef _WIN32
extern "C" {
    #include "compute_paths.h"
}
#else
#include "compute_paths.h"
#endif

namespace py = pybind11;

// Helper function to create a complex numpy array from two C arrays
py::array_t<std::complex<float>> make_complex_array(
    const float* real,
    const float* imag,
    size_t num_rx,
    size_t num_tx,
    size_t num_paths
) {
    size_t total_size = num_rx * num_tx * num_paths;
    std::complex<float>* complex_data = new std::complex<float>[total_size];

    for (size_t i = 0; i < total_size; ++i)
        complex_data[i] = std::complex<float>(real[i], imag[i]);

    auto capsule = py::capsule(complex_data, [](void* ptr) {
        delete[] static_cast<std::complex<float>*>(ptr);
    });
    
    return py::array_t<std::complex<float>>(
        {num_rx, num_tx, num_paths}, complex_data, capsule
    );
}

class PathsInfoPython {
public:
    size_t num_paths;
    py::array_t<float> directions_rx;
    py::array_t<float> directions_tx;
    py::array_t<std::complex<float>> a_te;
    py::array_t<std::complex<float>> a_tm;
    py::array_t<float> tau;
    py::array_t<float> freq_shift;

    PathsInfoPython(const PathsInfo& paths, size_t num_rx, size_t num_tx) {
        num_paths = paths.num_paths;

        auto capsule_deleter = [](void* ptr) { delete[] static_cast<float*>(ptr); };

        directions_rx = py::array_t<float>(
            std::vector<size_t>{num_rx, num_tx, (size_t)paths.num_paths, 3},
            paths.directions_rx,
            py::capsule(paths.directions_rx, capsule_deleter)
        );
        directions_tx = py::array_t<float>(
            std::vector<size_t>{num_rx, num_tx, (size_t)paths.num_paths, 3},
            paths.directions_tx,
            py::capsule(paths.directions_tx, capsule_deleter)
        );

        a_te = make_complex_array(paths.a_te_re, paths.a_te_im, num_rx, num_tx, paths.num_paths);
        a_tm = make_complex_array(paths.a_tm_re, paths.a_tm_im, num_rx, num_tx, paths.num_paths);

        tau = py::array_t<float>(
            {num_rx, num_tx, (size_t)paths.num_paths},
            paths.tau,
            py::capsule(paths.tau, capsule_deleter)
        );

        freq_shift = py::array_t<float>(
            {num_rx, num_tx, (size_t)paths.num_paths},
            paths.freq_shift,
            py::capsule(paths.freq_shift, capsule_deleter)
        );
    }
};

std::tuple<PathsInfoPython, PathsInfoPython>
compute_paths_wrapper(
    const std::string &mesh_filepath,
    py::array_t<float> rx_positions,
    py::array_t<float> tx_positions,
    py::array_t<float> rx_velocities,
    py::array_t<float> tx_velocities,
    float carrier_frequency,
    unsigned long int num_rx,
    unsigned long int num_tx,
    unsigned long int num_paths,
    unsigned long int num_bounces
) {
    // Prepare input arrays (this is a basic implementation, check shapes and memory layout)
    py::buffer_info rx_pos_info = rx_positions.request();
    py::buffer_info tx_pos_info = tx_positions.request();
    py::buffer_info rx_vel_info = rx_velocities.request();
    py::buffer_info tx_vel_info = tx_velocities.request();

    // Output
    PathsInfo los = {
        .num_paths = 1,
        .directions_rx = new float[num_rx * num_tx * 3],
        .directions_tx = new float[num_rx * num_tx * 3],
        .a_te_re = new float[num_rx * num_tx],
        .a_te_im = new float[num_rx * num_tx],
        .a_tm_re = new float[num_rx * num_tx],
        .a_tm_im = new float[num_rx * num_tx],
        .tau = new float[num_rx * num_tx],
        .freq_shift = new float[num_rx * num_tx]
    };
    PathsInfo scatter = {
        .num_paths = (uint32_t)(num_bounces * num_paths),
        .directions_rx = new float[num_rx * num_tx * num_bounces * num_paths * 3],
        .directions_tx = new float[num_rx * num_tx * num_bounces * num_paths * 3],
        .a_te_re = new float[num_rx * num_tx * num_bounces * num_paths],
        .a_te_im = new float[num_rx * num_tx * num_bounces * num_paths],
        .a_tm_re = new float[num_rx * num_tx * num_bounces * num_paths],
        .a_tm_im = new float[num_rx * num_tx * num_bounces * num_paths],
        .tau = new float[num_rx * num_tx * num_bounces * num_paths],
        .freq_shift = new float[num_rx * num_tx * num_bounces * num_paths]
    };

    // Call the C function
    compute_paths(
        mesh_filepath.c_str(),
        (const float*)rx_pos_info.ptr,  // Rx positions
        (const float*)tx_pos_info.ptr,  // Tx positions
        (const float*)rx_vel_info.ptr,  // Rx velocities
        (const float*)tx_vel_info.ptr,  // Tx velocities
        carrier_frequency,  // Carrier frequency in GHz
        (size_t)num_rx,
        (size_t)num_tx,
        (size_t)num_paths,
        (size_t)num_bounces,
        &los,
        &scatter
    );

    // Wrap the results into Python objects
    return std::make_tuple(
        PathsInfoPython(los, num_rx, num_tx),
        PathsInfoPython(scatter, num_rx, num_tx)
    );
}

PYBIND11_MODULE(rt, m) {
    py::class_<PathsInfoPython>(m, "PathsInfo")
        .def_readonly("num_paths", &PathsInfoPython::num_paths)
        .def_readonly("directions_rx", &PathsInfoPython::directions_rx)
        .def_readonly("directions_tx", &PathsInfoPython::directions_tx)
        .def_readonly("a_te", &PathsInfoPython::a_te)
        .def_readonly("a_tm", &PathsInfoPython::a_tm)
        .def_readonly("tau", &PathsInfoPython::tau)
        .def_readonly("freq_shift", &PathsInfoPython::freq_shift);

    // TODO write a proper docstring
    m.def("compute_paths", &compute_paths_wrapper, "Compute gains and delays",
          py::arg("mesh_filepath"),
          py::arg("rx_positions"),
          py::arg("tx_positions"),
          py::arg("rx_velocities"),
          py::arg("tx_velocities"),
          py::arg("carrier_frequency"),
          py::arg("num_rx"),
          py::arg("num_tx"),
          py::arg("num_paths"),
          py::arg("num_bounces"));
}

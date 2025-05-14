#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <complex>
#include <vector>
#include <iostream>

#ifdef _WIN32
extern "C" {
#endif
    #include "inc/compute_paths.h"
    #include "inc/scene.h"
    #include "inc/vec3.h"
    #include "inc/ray.h"
#ifdef _WIN32
}
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

class ChannelInfoPython {
public:
    size_t num_paths;
    py::array_t<float> directions_rx;
    py::array_t<float> directions_tx;
    py::array_t<std::complex<float>> a_te;
    py::array_t<std::complex<float>> a_tm;
    py::array_t<float> tau;
    py::array_t<float> freq_shift;

    ChannelInfoPython(const ChannelInfo& chanInfo, size_t num_rx, size_t num_tx) {
        num_paths = chanInfo.num_rays;

        auto capsule_deleter = [](void* ptr) { delete[] static_cast<float*>(ptr); };

        directions_rx = py::array_t<float>(
            std::vector<size_t>{num_rx, num_tx, (size_t)chanInfo.num_rays, 3},
            (float*)chanInfo.directions_rx,
            py::capsule(chanInfo.directions_rx, capsule_deleter)
        );
        directions_tx = py::array_t<float>(
            std::vector<size_t>{num_rx, num_tx, (size_t)chanInfo.num_rays, 3},
            (float*)chanInfo.directions_tx,
            py::capsule(chanInfo.directions_tx, capsule_deleter)
        );

        a_te = make_complex_array(
            chanInfo.a_te_re,
            chanInfo.a_te_im,
            num_rx,
            num_tx,
            chanInfo.num_rays
        );
        a_tm = make_complex_array(
            chanInfo.a_tm_re,
            chanInfo.a_tm_im,
            num_rx,
            num_tx,
            chanInfo.num_rays
        );

        tau = py::array_t<float>(
            {num_rx, num_tx, (size_t)chanInfo.num_rays},
            chanInfo.tau,
            py::capsule(chanInfo.tau, capsule_deleter)
        );

        freq_shift = py::array_t<float>(
            {num_rx, num_tx, (size_t)chanInfo.num_rays},
            chanInfo.freq_shift,
            py::capsule(chanInfo.freq_shift, capsule_deleter)
        );
    }
};

std::tuple<ChannelInfoPython, ChannelInfoPython>
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

    // Load the scene
    Scene scene = scene_load(mesh_filepath.c_str());

    // Output channel information
    ChannelInfo chanInfo_los = {
        .num_rays = 1,
        .directions_rx = new Vec3[num_rx * num_tx],
        .directions_tx = new Vec3[num_rx * num_tx],
        .a_te_re = new float[num_rx * num_tx],
        .a_te_im = new float[num_rx * num_tx],
        .a_tm_re = new float[num_rx * num_tx],
        .a_tm_im = new float[num_rx * num_tx],
        .tau = new float[num_rx * num_tx],
        .freq_shift = new float[num_rx * num_tx]
    };
    ChannelInfo chanInfo_scat = {
        .num_rays = (uint32_t)(num_bounces * num_paths),
        .directions_rx = new Vec3[num_rx * num_tx * num_bounces * num_paths],
        .directions_tx = new Vec3[num_rx * num_tx * num_bounces * num_paths],
        .a_te_re = new float[num_rx * num_tx * num_bounces * num_paths],
        .a_te_im = new float[num_rx * num_tx * num_bounces * num_paths],
        .a_tm_re = new float[num_rx * num_tx * num_bounces * num_paths],
        .a_tm_im = new float[num_rx * num_tx * num_bounces * num_paths],
        .tau = new float[num_rx * num_tx * num_bounces * num_paths],
        .freq_shift = new float[num_rx * num_tx * num_bounces * num_paths]
    };
    // Output rays information
    RaysInfo raysInfo_los = {
        .rays = new Ray[num_rx * num_tx],
        .rays_active = new uint8_t[num_rx * num_tx / 8 + 1]
    };
    RaysInfo raysInfo_scat = {
        .rays = new Ray[num_rx * num_tx * (num_bounces + 1) * num_paths],
        .rays_active = new uint8_t[num_rx * num_tx * (num_bounces + 1) * num_paths / 8 + 1]
    };

    // Call the C function
    compute_paths(
        &scene,  // Scene
        (Vec3*)rx_pos_info.ptr,  // Rx positions
        (Vec3*)tx_pos_info.ptr,  // Tx positions
        (Vec3*)rx_vel_info.ptr,  // Rx velocities
        (Vec3*)tx_vel_info.ptr,  // Tx velocities
        carrier_frequency,  // Carrier frequency in GHz
        (size_t)num_rx,
        (size_t)num_tx,
        (size_t)num_paths,
        (size_t)num_bounces,
        &chanInfo_los,  // Channel info for LOS
        &raysInfo_los,  // Rays info for LOS
        &chanInfo_scat,  // Channel info for scatter
        &raysInfo_scat  // Rays info for scatter
    );

    // Free the rays info (TODO remove when RaysInfo is made optional in compute_paths)
    delete[] raysInfo_los.rays;
    delete[] raysInfo_los.rays_active;
    delete[] raysInfo_scat.rays;
    delete[] raysInfo_scat.rays_active;

    // Free the scene
    free_scene(&scene);

    // Wrap the results into Python objects
    return std::make_tuple(
        ChannelInfoPython(chanInfo_los, num_rx, num_tx),
        ChannelInfoPython(chanInfo_scat, num_rx, num_tx)
    );
}

PYBIND11_MODULE(hermespy_rt, m) {
    py::class_<ChannelInfoPython>(m, "ChannelInfo")
        .def_readonly("num_paths", &ChannelInfoPython::num_paths)
        .def_readonly("directions_rx", &ChannelInfoPython::directions_rx)
        .def_readonly("directions_tx", &ChannelInfoPython::directions_tx)
        .def_readonly("a_te", &ChannelInfoPython::a_te)
        .def_readonly("a_tm", &ChannelInfoPython::a_tm)
        .def_readonly("tau", &ChannelInfoPython::tau)
        .def_readonly("freq_shift", &ChannelInfoPython::freq_shift);

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

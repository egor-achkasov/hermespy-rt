#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "compute_paths.h"

namespace py = pybind11;

std::tuple< py::array_t<float>, py::array_t<float> >
compute_paths_wrapper(
    const std::string &mesh_filepath,
    py::array_t<float> rx_positions,
    py::array_t<float> tx_positions,
    py::array_t<float> rx_velocities,
    py::array_t<float> tx_velocities,
    float carrier_frequency,
    int num_rx,
    int num_tx,
    int num_paths,
    int num_bounces
) {
    // Prepare input arrays (this is a basic implementation, check shapes and memory layout)
    py::buffer_info rx_pos_info = rx_positions.request();
    py::buffer_info tx_pos_info = tx_positions.request();
    py::buffer_info rx_vel_info = rx_velocities.request();
    py::buffer_info tx_vel_info = tx_velocities.request();

    // Output arrays
    // float *a_im = new float[num_paths];
    // float *a_re = new float[num_paths];
    float *a = new float[num_bounces * num_tx * num_paths * 3];
    float *tau = new float[num_paths];

    // Call the C function
    compute_paths(
        mesh_filepath.c_str(),
        (float*)rx_pos_info.ptr,  // Rx positions
        (float*)tx_pos_info.ptr,  // Tx positions
        (float*)rx_vel_info.ptr,  // Rx velocities
        (float*)tx_vel_info.ptr,  // Tx velocities
        carrier_frequency,
        (size_t)num_rx,
        (size_t)num_tx,
        (size_t)num_paths,
        (size_t)num_bounces,
        a,
        tau
    );

    // Convert output arrays into numpy arrays for easy use in Python
    //py::array_t<float> a_array = py::array_t<float>(num_paths, a);
    py::array_t<float> a_array = py::array_t<float>({num_bounces, num_tx, num_paths, 3}, a);
    py::array_t<float> tau_array = py::array_t<float>(num_paths, tau);

    // Deallocate arrays
    // delete[] a;
    // delete[] tau;

    // Return the results as a tuple (gains, delays)
    return std::make_tuple(a_array, tau_array);
}

PYBIND11_MODULE(rt, m) {
    m.def("compute_paths", &compute_paths_wrapper, "Compute gains and delays in PLY scene",
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

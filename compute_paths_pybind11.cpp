#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#ifdef _WIN32
extern "C" {
    #include "compute_paths.h"
}
#else
#include "compute_paths.h"
#endif

namespace py = pybind11;

std::tuple<
  py::array_t<float>, py::array_t<float>,  // a_te_re, a_te_im
  py::array_t<float>, py::array_t<float>,  // a_tm_re, a_tm_im
  py::array_t<float> >                     // tau
compute_paths_wrapper(
    const std::string &mesh_filepath,
    py::array_t<float> rx_positions,
    py::array_t<float> tx_positions,
    py::array_t<float> rx_velocities,
    py::array_t<float> tx_velocities,
    float carrier_frequency,
    float sampling_frequency,
    int num_rx,
    int num_tx,
    int num_paths,
    int num_bounces,
    int num_samples
) {
    // Prepare input arrays (this is a basic implementation, check shapes and memory layout)
    py::buffer_info rx_pos_info = rx_positions.request();
    py::buffer_info tx_pos_info = tx_positions.request();
    py::buffer_info rx_vel_info = rx_velocities.request();
    py::buffer_info tx_vel_info = tx_velocities.request();

    // Output
    size_t num_paths_out;
    float *a_te_re, *a_te_im, *a_tm_re, *a_tm_im, *tau;

    // Call the C function
    compute_paths(
        mesh_filepath.c_str(),
        (const float*)rx_pos_info.ptr,  // Rx positions
        (const float*)tx_pos_info.ptr,  // Tx positions
        (const float*)rx_vel_info.ptr,  // Rx velocities
        (const float*)tx_vel_info.ptr,  // Tx velocities
        carrier_frequency,  // Carrier frequency in GHz
        sampling_frequency, // Sampling frequency in Hz
        (size_t)num_rx,
        (size_t)num_tx,
        (size_t)num_paths,
        (size_t)num_bounces,
        (size_t)num_samples,
        &num_paths_out,
        a_te_re, a_te_im, a_tm_re, a_tm_im,  // Gains
        tau // Delays
    );

    // Convert output arrays into numpy arrays for easy use in Python
    // TODO prevent copying data
    // TODO construct np.complex64 array for a
    py::array_t<float> a_te_re_array = py::array_t<float>({num_rx, num_tx, (int)num_paths_out}, a_te_re);
    py::array_t<float> a_te_im_array = py::array_t<float>({num_rx, num_tx, (int)num_paths_out}, a_te_im);
    py::array_t<float> a_tm_re_array = py::array_t<float>({num_rx, num_tx, (int)num_paths_out}, a_tm_re);
    py::array_t<float> a_tm_im_array = py::array_t<float>({num_rx, num_tx, (int)num_paths_out}, a_tm_im);
    py::array_t<float> tau_array = py::array_t<float>({num_rx, num_tx, (int)num_paths_out}, tau);

    // Deallocate arrays
    delete[] a_te_re;
    delete[] a_te_im;
    delete[] a_tm_re;
    delete[] a_tm_im;
    delete[] tau;

    // Return the results as a tuple (gains, delays)
    return std::make_tuple(
      a_te_re_array, a_te_im_array,
      a_tm_re_array, a_tm_im_array,
      tau_array);
}

PYBIND11_MODULE(rt, m) {
    m.def("compute_paths", &compute_paths_wrapper, "Compute gains and delays in PLY scene",
          py::arg("mesh_filepath"),
          py::arg("rx_positions"),
          py::arg("tx_positions"),
          py::arg("rx_velocities"),
          py::arg("tx_velocities"),
          py::arg("carrier_frequency"),
          py::arg("sampling_frequency"),
          py::arg("num_rx"),
          py::arg("num_tx"),
          py::arg("num_paths"),
          py::arg("num_bounces"),
          py::arg("num_samples"));

    // Scene filepaths
    m.def("get_scene_fp_box",
        []() { return __FILE__ + std::string("/scenes/box.ply"); }
    );
}

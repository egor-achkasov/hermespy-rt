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

// Helper function to create a numpy array from a C array
py::array_t<float> make_py_array(std::vector<size_t> shape, float *data) {
    return py::array_t<float>(
        //shape, data, py::capsule(data, [](float *ptr) { delete[] ptr; })
        shape, data
    );
}

std::tuple<
  py::array_t<float>,                      // directions_los
  py::array_t<float>, py::array_t<float>,  // a_te_re_los, a_te_im_los
  py::array_t<float>, py::array_t<float>,  // a_tm_re_los, a_tm_im_los
  py::array_t<float>,                      // tau_los
  py::array_t<float>,                      // directions_scat
  py::array_t<float>, py::array_t<float>,  // a_te_re_scat, a_te_im_scat
  py::array_t<float>, py::array_t<float>,  // a_tm_re_scat, a_tm_im_scat
  py::array_t<float>                       // tau_scat
>
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

    // Output
    float *directions_los = new float[num_rx * num_tx * 3];
    float *a_te_re_los = new float[num_rx * num_tx];
    float *a_te_im_los = new float[num_rx * num_tx];
    float *a_tm_re_los = new float[num_rx * num_tx];
    float *a_tm_im_los = new float[num_rx * num_tx];
    float *tau_los = new float[num_rx * num_tx];
    float *directions_scat = new float[num_rx * num_tx * num_paths * 3];
    float *a_te_re_scat = new float[num_bounces * num_rx * num_tx * num_paths];
    float *a_te_im_scat = new float[num_bounces * num_rx * num_tx * num_paths];
    float *a_tm_re_scat = new float[num_bounces * num_rx * num_tx * num_paths];
    float *a_tm_im_scat = new float[num_bounces * num_rx * num_tx * num_paths];
    float *tau_scat = new float[num_bounces * num_rx * num_tx * num_paths];

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
        // LoS outputs
        directions_los,
        a_te_re_los,
        a_te_im_los,
        a_tm_re_los,
        a_tm_im_los,
        tau_los,
        // Scatter outputs
        directions_scat,
        a_te_re_scat,
        a_te_im_scat,
        a_tm_re_scat,
        a_tm_im_scat,
        tau_scat
    );

    // Cast to numpy arrays and return
    return std::make_tuple(
      make_py_array({num_rx, num_tx, 3}, directions_los),
      make_py_array({num_rx, num_tx}, a_te_re_los),
      make_py_array({num_rx, num_tx}, a_te_im_los),
      make_py_array({num_rx, num_tx}, a_tm_re_los),
      make_py_array({num_rx, num_tx}, a_tm_im_los),
      make_py_array({num_rx, num_tx}, tau_los),
      make_py_array({num_rx, num_tx, num_paths, 3}, directions_scat),
      make_py_array({num_bounces, num_rx, num_tx, num_paths}, a_te_re_scat),
      make_py_array({num_bounces, num_rx, num_tx, num_paths}, a_te_im_scat),
      make_py_array({num_bounces, num_rx, num_tx, num_paths}, a_tm_re_scat),
      make_py_array({num_bounces, num_rx, num_tx, num_paths}, a_tm_im_scat),
      make_py_array({num_bounces, num_rx, num_tx, num_paths}, tau_scat)
    );
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

    // Scene filepaths
    // TODO
    m.def("get_scene_fp_box",
        []() { return __FILE__ + std::string("/scenes/box.ply"); }
    );
}

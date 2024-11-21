/* This file contains a python binding for compute_paths function. */

#include <Python.h>
#include <numpy/arrayobject.h>  // Numpy C API
#include "compute_paths.h"      // Assuming this header contains the compute_paths function declaration

// Wrapper function for Python to call compute_paths
static PyObject* py_compute_paths(PyObject* self, PyObject* args) {
    const char* scene_filepath;
    PyObject* rx_positions_obj;
    PyObject* tx_positions_obj;
    PyObject* rx_velocities_obj;
    PyObject* tx_velocities_obj;
    float carrier_frequency;
    int num_rx, num_tx, num_paths, num_bounces;
    
    if (!PyArg_ParseTuple(args, "sO&O&O&O&fiiii",
                          &scene_filepath,
                          &rx_positions_obj,
                          &tx_positions_obj,
                          &rx_velocities_obj,
                          &tx_velocities_obj,
                          &carrier_frequency,
                          &num_rx,
                          &num_tx,
                          &num_paths,
                          &num_bounces))
        return NULL;

    if (!PyArray_Check(rx_positions_obj) || !PyArray_Check(tx_positions_obj) ||
        !PyArray_Check(rx_velocities_obj) || !PyArray_Check(tx_velocities_obj)) {
        PyErr_SetString(PyExc_TypeError, "All inputs must be numpy arrays.");
        return NULL;
    }

    npy_float32* rx_positions = (npy_float32*)PyArray_DATA((PyArrayObject*)rx_positions_obj);
    npy_float32* tx_positions = (npy_float32*)PyArray_DATA((PyArrayObject*)tx_positions_obj);
    npy_float32* rx_velocities = (npy_float32*)PyArray_DATA((PyArrayObject*)rx_velocities_obj);
    npy_float32* tx_velocities = (npy_float32*)PyArray_DATA((PyArrayObject*)tx_velocities_obj);

    // Output arrays for gains and delays
    float* a = (float*)malloc(num_paths * sizeof(float));
    int32_t* tau = (int32_t*)malloc(num_paths * sizeof(int32_t));

    if (!a || !tau) {
        PyErr_NoMemory();
        return NULL;
    }

    // Call the C compute_paths function
    compute_paths(scene_filepath, 
                  rx_positions, tx_positions, 
                  rx_velocities, tx_velocities, 
                  carrier_frequency,
                  num_rx, num_tx, num_paths, num_bounces,
                  a, tau);

    // Prepare the output tuple (gains, delays)
    //npy_intp dims[1] = {num_paths};
    npy_intp dims[] = {num_bounces, num_tx, num_paths};
    PyObject* gains_array = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT32, a);
    PyObject* delays_array = PyArray_SimpleNewFromData(1, dims, NPY_INT32, tau);

    // Return the tuple of numpy arrays
    PyObject* result = PyTuple_Pack(2, gains_array, delays_array);

    // Clean up
    Py_DECREF(gains_array);
    Py_DECREF(delays_array);

    return result;
}

// Method table
/* TODO add args and returns */
static PyMethodDef module_methods[] = {
    {"compute_paths", py_compute_paths, METH_VARARGS, "Compute the paths between transmitters and receivers."},
    {NULL, NULL, 0, NULL}
};

// Module definition
static struct PyModuleDef module_definition = {
    PyModuleDef_HEAD_INIT,
    "compute_paths_module",  // Module name
    "C extension for computing paths.",  // Module docstring
    -1,
    module_methods  // Method table
};

// Module initialization function
PyMODINIT_FUNC PyInit_compute_paths_module(void) {
    import_array();  // Initialize the numpy API
    return PyModule_Create(&module_definition);
}

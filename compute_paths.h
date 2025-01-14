#ifndef COMPUTE_PATHS_H
#define COMPUTE_PATHS_H

#include <stddef.h> /* for size_t */

#define IN
#define OUT

/** Compute gains and delays between in a Mitsuba scene.
 * 
 * Scene must be defined in a specific PLY format. See README for details.
 * The meshes must be present in the same directory as the .xml file,
 * under meshes/ directory in .ply files in PLY format.
 * Meshes filenames must be the same as the name of the object in the .xml file.
 * 
 * \param mesh_filepath path to a mesh .ply file
 * \param rx_positions receiver positions, shape (num_rx, 3)
 * \param tx_positions transmitter positions, shape (num_tx, 3)
 * \param rx_velocities receiver velocities, shape (num_rx, 3)
 * \param tx_velocities transmitter velocities, shape (num_tx, 3)
 * \param carrier_frequency carrier frequency in GHz. Must be > 0.0
 * \param num_rx number of receivers. Must be > 0
 * \param num_tx number of transmitters. Must be > 0
 * \param num_paths number of paths to compute. Must be > 0
 * \param num_bounces number of bounces to compute. Must be > 0
 * \param num_samples number of samples in the original signal. Must be > 0
 * \param a_te_re output array of real parts of transverse electric gains, shape (num_rx, num_tx, num_paths)
 * \param a_te_im output array of imaginary parts of transverse electric gains, shape (num_rx, num_tx, num_paths)
 * \param a_tm_re output array of real parts of transverse magnetic gains, shape (num_rx, num_tx, num_paths)
 * \param a_tm_im output array of imaginary parts of transverse magnetic gains, shape (num_rx, num_tx, num_paths)
 * \param tau output array of delays in seconds, shape (num_rx, num_tx, num_paths)
*/
void compute_paths(
    IN const char *mesh_filepath,  /* path to the mesh file */
    IN const float *rx_positions,  /* shape (num_rx, 3) */
    IN const float *tx_positions,  /* shape (num_tx, 3) */
    IN const float *rx_velocities, /* shape (num_rx, 3) */
    IN const float *tx_velocities, /* shape (num_tx, 3) */
    IN float carrier_frequency,    /* > 0.0 (IN GHz!) */
    IN size_t num_rx,              /* number of receivers */
    IN size_t num_tx,              /* number of transmitters */
    IN size_t num_paths,           /* number of paths */
    IN size_t num_bounces,         /* number of bounces */
    IN size_t num_samples,         /* number of samples */
    OUT float *a_te_re,            /* output array real parts of TE gains (num_rx, num_tx, num_paths) */
    OUT float *a_te_im,            /* output array imaginary parts of TE gains (num_rx, num_tx, num_paths) */
    OUT float *a_tm_re,            /* output array real parts of TM gains (num_rx, num_tx, num_paths) */
    OUT float *a_tm_im,            /* output array imaginary parts of TM gains (num_rx, num_tx, num_paths) */
    OUT float *tau                 /* output array of delays (num_rx, num_tx, num_paths) */
);

#endif  /* COMPUTE_PATHS_H */

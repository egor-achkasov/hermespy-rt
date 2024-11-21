#ifndef COMPUTE_PATHS_H
#define COMPUTE_PATHS_H

#include <stdint.h>  /* for int32_t */

#include "common.h"

/** Compute gains and delays between in a Mitsuba scene.
 * 
 * Scene must be defined in a Mitsuba .xml JSON format.
 * The meshes must be present in the same directory as the .xml file,
 * under meshes/ directory in .ply files in PLY format.
 * Meshes filenames must be the same as the name of the object in the .xml file.
 * 
 * \param scene_filepath path to a Mitsuba scene .xml file
 * \param rx_positions receiver positions, shape (num_rx, 3)
 * \param tx_positions transmitter positions, shape (num_tx, 3)
 * \param rx_velocities receiver velocities, shape (num_rx, 3)
 * \param tx_velocities transmitter velocities, shape (num_tx, 3)
 * \param carrier_frequency carrier frequency in Hz. Must be > 0.0
 * \param num_rx number of receivers. Must be > 0
 * \param num_tx number of transmitters. Must be > 0
 * \param num_paths number of paths to compute. Must be > 0
 * \param a output array of gains, shape (num_paths,)
 * \param tau output array of delays, shape (num_paths,)
*/
void compute_paths(
    IN const char *scene_filepath, /* path to the scene file */
    IN const float *rx_positions,  /* shape (num_rx, 3) */
    IN const float *tx_positions,  /* shape (num_tx, 3) */
    IN const float *rx_velocities, /* shape (num_rx, 3) */
    IN const float *tx_velocities, /* shape (num_tx, 3) */
    IN float carrier_frequency,    /* > 0.0 */
    IN int num_rx,                 /* number of receivers */
    IN int num_tx,                 /* number of transmitters */
    IN int num_paths,              /* number of paths */
    IN int num_bounces,            /* number of bounces */
    OUT float *a,                  /* output array of gains (num_paths,) */
    OUT int32_t *tau               /* output array of delays (num_paths,) */
);

#endif  /* COMPUTE_PATHS_H */

#ifndef COMPUTE_PATHS_H
#define COMPUTE_PATHS_H

#include "scene.h" /* for Scene */

#include <stddef.h> /* for size_t */
#include <stdint.h> /* for uint32_t */

#define IN
#define OUT

/** Raytracing information returned by compute_paths for each path type (LoS, scatter) */
typedef struct {
    uint32_t num_paths;
    float *directions_rx; /* shape (num_rx, num_tx, num_paths, 3) */
    float *directions_tx; /* shape (num_rx, num_tx, num_paths, 3) */
    float *a_te_re; /* shape (num_rx, num_tx, num_paths) */
    float *a_te_im; /* shape (num_rx, num_tx, num_paths) */
    float *a_tm_re; /* shape (num_rx, num_tx, num_paths) */
    float *a_tm_im; /* shape (num_rx, num_tx, num_paths) */
    float *tau; /* shape (num_rx, num_tx, num_paths) */
    float *freq_shift; /* Hz, shape (num_rx, num_tx, num_paths) */
} PathsInfo;

/** Compute gains and delays between tx and rx in a 3D scene.
 * 
 * Scene must be defined in a specific format.
 * The format is defined in scene.h.
 * The HRT is a binary file cointaining a dereferenced Scene structure.
 * See README for details.
 * 
 * The output PathInfo structure is defined in compute_paths.h
 * 
 * The output parameters are allocated by the caller (including the fields).
 * 
 * \param scene pointer to a loaded scene in HRT format (see scene.h)
 * \param rx_pos receiver positions, shape (num_rx, 3)
 * \param tx_pos transmitter positions, shape (num_tx, 3)
 * \param rx_vel receiver velocities, shape (num_rx, 3)
 * \param tx_vel transmitter velocities, shape (num_tx, 3)
 * \param carrier_frequency carrier frequency in GHz. Must be > 0.0
 * \param num_rx number of receivers. Must be > 0
 * \param num_tx number of transmitters. Must be > 0
 * \param num_paths number of paths to compute. Must be > 0
 * \param num_bounces number of bounces to compute. Must be > 0
 * 
 * \param los output Line-of-Sight results. los->num_paths = 1
 * \param scatter output scatter results. scatter->num_paths = num_bounces * num_paths
*/
void compute_paths(
    IN Scene *scene,            /* Pointer to a loaded scene */
    IN Vec3 *rx_pos,            /* shape (num_rx, 3) */
    IN Vec3 *tx_pos,            /* shape (num_tx, 3) */
    IN Vec3 *rx_vel,            /* shape (num_rx, 3) */
    IN Vec3 *tx_vel,            /* shape (num_tx, 3) */
    IN float carrier_frequency, /* > 0.0 (IN GHz!) */
    IN size_t num_rx,           /* number of receivers */
    IN size_t num_tx,           /* number of transmitters */
    IN size_t num_paths,        /* number of paths */
    IN size_t num_bounces,      /* number of bounces */
    OUT PathsInfo *los,         /* output LoS information */
    OUT PathsInfo *scatter      /* output scatter information */
);

#endif /* COMPUTE_PATHS_H */

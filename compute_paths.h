#ifndef COMPUTE_PATHS_H
#define COMPUTE_PATHS_H

#include <stddef.h> /* for size_t */

#define IN
#define OUT

/** Compute gains and delays between in a Mitsuba scene.
 * 
 * Scene must be defined in a specific PLY format. See README for details.
 * 
 * \param mesh_filepath path to a mesh .ply file
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
 * Outputs:
 * 
 * LoS:
 * \param directions_los output array of directions of incidence for LoS, shape (num_rx, num_tx, 3)
 * \param a_te_re_los output array of real parts of transverse electric gains for LoS, shape (num_rx, num_tx)
 * \param a_te_im_los output array of imaginary parts of transverse electric gains for LoS, shape (num_rx, num_tx)
 * \param a_tm_re_los output array of real parts of transverse magnetic gains for LoS, shape (num_rx, num_tx)
 * \param a_tm_im_los output array of imaginary parts of transverse magnetic gains for LoS, shape (num_rx, num_tx)
 * \param tau_los output array of delays for LoS in seconds, shape (num_rx, num_tx)
 * 
 * Scatter:
 * \param directions_scat output array of directions of incidence for scatter, shape (num_rx, num_tx, num_paths, 3)
 * \param a_te_re_scat output array of real parts of transverse electric gains for scatter, shape (num_bounces, num_rx, num_tx, num_paths)
 * \param a_te_im_scat output array of imaginary parts of transverse electric gains for scatter, shape (num_bounces, num_rx, num_tx, num_paths)
 * \param a_tm_re_scat output array of real parts of transverse magnetic gains for scatter, shape (num_bounces, num_rx, num_tx, num_paths)
 * \param a_tm_im_scat output array of imaginary parts of transverse magnetic gains for scatter, shape (num_bounces, num_rx, num_tx, num_paths)
 * \param tau_scat output array of delays for scatter in seconds, shape (num_bounces, num_rx, num_tx, num_paths)
*/
void compute_paths(
    IN const char *mesh_filepath,   /* path to the mesh file */
    IN const float *rx_pos,         /* shape (num_rx, 3) */
    IN const float *tx_pos,         /* shape (num_tx, 3) */
    IN const float *rx_vel,         /* shape (num_rx, 3) */
    IN const float *tx_vel,         /* shape (num_tx, 3) */
    IN float carrier_frequency,     /* > 0.0 (IN GHz!) */
    IN size_t num_rx,               /* number of receivers */
    IN size_t num_tx,               /* number of transmitters */
    IN size_t num_paths,            /* number of paths */
    IN size_t num_bounces,          /* number of bounces */
    /* LoS */
    OUT float *directions_los,      /* output array of directions of incidence (num_rx, num_tx, 3) */
    OUT float *a_te_re_los,         /* output array real parts of TE gains (num_rx, num_tx) */
    OUT float *a_te_im_los,         /* output array imaginary parts of TE gains (num_rx, num_tx) */
    OUT float *a_tm_re_los,         /* output array real parts of TM gains (num_rx, num_tx) */
    OUT float *a_tm_im_los,         /* output array imaginary parts of TM gains (num_rx, num_tx) */
    OUT float *tau_los,             /* output array of delays (num_rx, num_tx) */
    /* Scatter */
    OUT float *directions_scat,     /* output array of directions of incidence (num_rx, num_tx, num_paths, 3) */
    OUT float *a_te_re_scat,        /* output array real parts of TE gains (num_bounces, num_rx, num_tx, num_paths) */
    OUT float *a_te_im_scat,        /* output array imaginary parts of TE gains (num_bounces, num_rx, num_tx, num_paths) */
    OUT float *a_tm_re_scat,        /* output array real parts of TM gains (num_bounces, num_rx, num_tx, num_paths) */
    OUT float *a_tm_im_scat,        /* output array imaginary parts of TM gains (num_bounces, num_rx, num_tx, num_paths) */
    OUT float *tau_scat             /* output array of delays (num_bounces, num_rx, num_tx, num_paths) */
);

#endif /* COMPUTE_PATHS_H */

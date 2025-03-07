#ifndef COMPUTE_PATHS_H
#define COMPUTE_PATHS_H

#include <stddef.h> /* for size_t */

#define IN
#define OUT

/** Compute gains and delays between tx and rx in a 3D scene.
 * 
 * Scene must be defined in a specific format.
 * The format is defined in scene.h.
 * The HRT is a binary file cointaining a dereferenced Scene structure.
 * See README for details.
 * 
 * \param scene_filepath path to a scene .hrt file
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
 * \param directions_rx_scat output array of directions of incidence for rx for scatter, shape (num_rx, num_tx, num_bounces * num_paths, 3)
 * \param directions_tx_scat output array of directions of incidence for tx for scatter, shape (num_paths, 3)
 * \param a_te_re_scat output array of real parts of transverse electric gains for scatter, shape (num_rx, num_tx, num_bounces * num_paths)
 * \param a_te_im_scat output array of imaginary parts of transverse electric gains for scatter, shape (num_rx, num_tx, num_bounces * num_paths)
 * \param a_tm_re_scat output array of real parts of transverse magnetic gains for scatter, shape (num_rx, num_tx, num_bounces * num_paths)
 * \param a_tm_im_scat output array of imaginary parts of transverse magnetic gains for scatter, shape (num_rx, num_tx, num_bounces * num_paths)
 * \param tau_scat output array of delays for scatter in seconds, shape (num_rx, num_tx, num_bounces * num_paths)
*/
void compute_paths(
    IN const char *scene_filepath,   /* path to the scene file */
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
    OUT float *directions_los,   /* output array of directions of incidence for rx (num_rx, num_tx, 3) */
    OUT float *a_te_re_los,         /* output array real parts of TE gains (num_rx, num_tx) */
    OUT float *a_te_im_los,         /* output array imaginary parts of TE gains (num_rx, num_tx) */
    OUT float *a_tm_re_los,         /* output array real parts of TM gains (num_rx, num_tx) */
    OUT float *a_tm_im_los,         /* output array imaginary parts of TM gains (num_rx, num_tx) */
    OUT float *tau_los,             /* output array of delays (num_rx, num_tx) */
    /* Scatter */
    OUT float *directions_rx_scat,  /* output array of directions of incidence for rx (num_rx, num_tx, num_bounces * num_paths, 3) */
    OUT float *directions_tx_scat,  /* output array of directions of incidence for tx (num_paths, 3) */
    OUT float *a_te_re_scat,        /* output array real parts of TE gains (num_rx, num_tx, num_ounces * num_paths) */
    OUT float *a_te_im_scat,        /* output array imaginary parts of TE gains (num_rx, num_tx, num_ounces * num_paths) */
    OUT float *a_tm_re_scat,        /* output array real parts of TM gains (num_rx, num_tx, num_ounces * num_paths) */
    OUT float *a_tm_im_scat,        /* output array imaginary parts of TM gains (num_rx, num_tx, num_ounces * num_paths) */
    OUT float *tau_scat             /* output array of delays (num_rx, num_tx, num_ounces * num_paths) */
);

#endif /* COMPUTE_PATHS_H */

#include "compute_paths.h"

#include <stdio.h> /* for fprintf */

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <path_to_mesh.ply>\n", argv[0]);
        return 1;
    }

    float rx_positions[3] = {0.0, 0.0, 2.5};
    float tx_positions[3] = {0.0, 0.0, 2.5};
    float rx_velocities[3] = {0.0, 0.0, 0.0};
    float tx_velocities[3] = {0.0, 0.0, 0.0};
    float carrier_frequency = 2.4e9; /* 2.4 GHz */
    size_t num_paths = 10000;
    size_t num_bounces = 3;
    float *a = NULL;
    float *tau = NULL;
    compute_paths(
        argv[1],
        rx_positions,
        tx_positions,
        rx_velocities,
        tx_velocities,
        carrier_frequency,
        1, 1, num_paths, num_bounces,
        a, tau);
}

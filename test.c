#include "compute_paths.h"

#include <stdio.h> /* for fprintf */

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <path_to_scene.xml>\n", argv[0]);
        return 1;
    }

    float rx_positions[3] = {0.0, 0.0, 0.0};
    float tx_positions[3] = {0.0, 0.0, 0.0};
    float rx_velocities[3] = {0.0, 0.0, 0.0};
    float tx_velocities[3] = {0.0, 0.0, 0.0};
    float carrier_frequency = 2.4e9; /* 2.4 GHz */
    int num_paths = 10000;
    float *a = NULL;
    int32_t *tau = NULL;
    compute_paths(
        argv[1],
        rx_positions,
        tx_positions,
        rx_velocities,
        tx_velocities,
        carrier_frequency,
        1, 1, num_paths, 3,
        a, tau);
}

#include "compute_paths.h"

#include <stdio.h> /* for fprintf */
#include <stdlib.h> /* for malloc */

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
  float carrier_frequency = 3.0; /* 3 GHz */
  size_t num_rx = 1;
  size_t num_tx = 1;
  size_t num_paths = 10000;
  size_t num_bounces = 3;
  size_t num_samples = 150;
  float *a_te_re = (float*)malloc(num_rx * num_tx * num_paths * sizeof(float));
  float *a_te_im = (float*)malloc(num_rx * num_tx * num_paths * sizeof(float));
  float *a_tm_re = (float*)malloc(num_rx * num_tx * num_paths * sizeof(float));
  float *a_tm_im = (float*)malloc(num_rx * num_tx * num_paths * sizeof(float));
  float *tau = (float*)malloc(num_rx * num_tx * num_paths * sizeof(float));
  compute_paths(
    argv[1],
    rx_positions,
    tx_positions,
    rx_velocities,
    tx_velocities,
    carrier_frequency,
    num_rx, num_tx, num_paths, num_bounces, num_samples,
    a_te_re, a_te_im, a_tm_re, a_tm_im,
    tau);

  /* Save the results */
  FILE *f = fopen("a_te_re.bin", "wb");
  fwrite(a_te_re, sizeof(float), num_rx * num_tx * num_paths, f);
  fclose(f);
  f = fopen("a_te_im.bin", "wb");
  fwrite(a_te_im, sizeof(float), num_rx * num_tx * num_paths, f);
  fclose(f);
  f = fopen("a_tm_re.bin", "wb");
  fwrite(a_tm_re, sizeof(float), num_rx * num_tx * num_paths, f);
  fclose(f);
  f = fopen("a_tm_im.bin", "wb");
  fwrite(a_tm_im, sizeof(float), num_rx * num_tx * num_paths, f);
  fclose(f);
  f = fopen("tau.bin", "wb");
  fwrite(tau, sizeof(float), num_rx * num_tx * num_paths, f);
  fclose(f);
}

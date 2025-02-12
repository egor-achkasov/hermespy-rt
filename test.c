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
  float *directions_los = (float*)malloc(num_rx * num_tx * 3 * sizeof(float));
  float *a_te_re_los = (float*)malloc(num_rx * num_tx * sizeof(float));
  float *a_te_im_los = (float*)malloc(num_rx * num_tx * sizeof(float));
  float *a_tm_re_los = (float*)malloc(num_rx * num_tx * sizeof(float));
  float *a_tm_im_los = (float*)malloc(num_rx * num_tx * sizeof(float));
  float *tau_los = (float*)malloc(num_rx * num_tx * sizeof(float));
  float *directions_scat = (float*)malloc(num_rx * num_tx * num_paths * 3 * sizeof(float));
  float *a_te_re_scat = (float*)malloc(num_rx * num_tx * num_bounces * num_paths * sizeof(float));
  float *a_te_im_scat = (float*)malloc(num_rx * num_tx * num_bounces * num_paths * sizeof(float));
  float *a_tm_re_scat = (float*)malloc(num_rx * num_tx * num_bounces * num_paths * sizeof(float));
  float *a_tm_im_scat = (float*)malloc(num_rx * num_tx * num_bounces * num_paths * sizeof(float));
  float *tau_scat = (float*)malloc(num_rx * num_tx * num_bounces * num_paths * sizeof(float));
  compute_paths(
    argv[1],
    rx_positions,
    tx_positions,
    rx_velocities,
    tx_velocities,
    carrier_frequency,
    num_rx, num_tx, num_paths, num_bounces,
    directions_los, a_te_re_los, a_te_im_los, a_tm_re_los, a_tm_im_los, tau_los,
    directions_scat, a_te_re_scat, a_te_im_scat, a_tm_re_scat, a_tm_im_scat, tau_scat
  );

  /* Save the results */

  FILE *f;
  #define WRITE_BIN(name, data, size) \
    f = fopen(name, "wb"); \
    fwrite(data, sizeof(float), size, f); \
    fclose(f);

  WRITE_BIN("directions_los.bin", directions_los, num_rx * num_tx * 3);
  WRITE_BIN("a_te_re_los.bin", a_te_re_los, num_rx * num_tx);
  WRITE_BIN("a_te_im_los.bin", a_te_im_los, num_rx * num_tx);
  WRITE_BIN("a_tm_re_los.bin", a_tm_re_los, num_rx * num_tx);
  WRITE_BIN("a_tm_im_los.bin", a_tm_im_los, num_rx * num_tx);
  WRITE_BIN("tau_los.bin", tau_los, num_rx * num_tx);
  WRITE_BIN("directions_scat.bin", directions_scat, num_rx * num_tx * num_paths * 3);
  WRITE_BIN("a_te_re_scat.bin", a_te_re_scat, num_rx * num_tx * num_bounces * num_paths);
  WRITE_BIN("a_te_im_scat.bin", a_te_im_scat, num_rx * num_tx * num_bounces * num_paths);
  WRITE_BIN("a_tm_re_scat.bin", a_tm_re_scat, num_rx * num_tx * num_bounces * num_paths);
  WRITE_BIN("a_tm_im_scat.bin", a_tm_im_scat, num_rx * num_tx * num_bounces * num_paths);
  WRITE_BIN("tau_scat.bin", tau_scat, num_rx * num_tx * num_bounces * num_paths);
}

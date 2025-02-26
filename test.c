#include "compute_paths.h"
#include "test.h"

#include <stdio.h> /* for fprintf */
#include <stdlib.h> /* for malloc */

int main(int argc, char **argv)
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <path_to_mesh.ply>\n", argv[0]);
    return 1;
  }

  float rx_positions[3] = {0.0, 0.0, .5};
  float tx_positions[3] = {0.0, 0.0, .5};
  float rx_velocities[3] = {0.0, 0.0, 0.0};
  float tx_velocities[3] = {0.0, 0.0, 0.0};
  float carrier_frequency = 3.0; /* 3 GHz */
  float *directions_los = (float*)malloc(g_numRx * g_numTx * 3 * sizeof(float));
  float *a_te_re_los = (float*)malloc(g_numRx * g_numTx * sizeof(float));
  float *a_te_im_los = (float*)malloc(g_numRx * g_numTx * sizeof(float));
  float *a_tm_re_los = (float*)malloc(g_numRx * g_numTx * sizeof(float));
  float *a_tm_im_los = (float*)malloc(g_numRx * g_numTx * sizeof(float));
  float *tau_los = (float*)malloc(g_numRx * g_numTx * sizeof(float));
  float *directions_rx_scat = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * 3 * sizeof(float));
  float *directions_tx_scat = (float*)malloc(g_numPaths * 3 * sizeof(float));
  float *a_te_re_scat = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float));
  float *a_te_im_scat = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float));
  float *a_tm_re_scat = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float));
  float *a_tm_im_scat = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float));
  float *tau_scat = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float));
  compute_paths(
    argv[1],
    rx_positions,
    tx_positions,
    rx_velocities,
    tx_velocities,
    carrier_frequency,
    g_numRx, g_numTx, g_numPaths, g_numBounces,
    directions_los, a_te_re_los, a_te_im_los, a_tm_re_los, a_tm_im_los, tau_los,
    directions_rx_scat, directions_rx_scat,
    a_te_re_scat, a_te_im_scat, a_tm_re_scat, a_tm_im_scat, tau_scat
  );

  /* Save the results */

  FILE *f;
  #define WRITE_BIN(name, data, size) \
    f = fopen(name, "wb"); \
    fwrite(data, sizeof(float), size, f); \
    fclose(f);

  WRITE_BIN("directions_los.bin", directions_los, g_numRx * g_numTx * 3);
  WRITE_BIN("a_te_re_los.bin", a_te_re_los, g_numRx * g_numTx);
  WRITE_BIN("a_te_im_los.bin", a_te_im_los, g_numRx * g_numTx);
  WRITE_BIN("a_tm_re_los.bin", a_tm_re_los, g_numRx * g_numTx);
  WRITE_BIN("a_tm_im_los.bin", a_tm_im_los, g_numRx * g_numTx);
  WRITE_BIN("tau_los.bin", tau_los, g_numRx * g_numTx);
  WRITE_BIN("directions_rx_scat.bin", directions_rx_scat, g_numRx * g_numTx * g_numBounces * g_numPaths * 3);
  WRITE_BIN("directions_tx_scat.bin", directions_tx_scat, g_numPaths * 3);
  WRITE_BIN("a_te_re_scat.bin", a_te_re_scat, g_numRx * g_numTx * g_numBounces * g_numPaths);
  WRITE_BIN("a_te_im_scat.bin", a_te_im_scat, g_numRx * g_numTx * g_numBounces * g_numPaths);
  WRITE_BIN("a_tm_re_scat.bin", a_tm_re_scat, g_numRx * g_numTx * g_numBounces * g_numPaths);
  WRITE_BIN("a_tm_im_scat.bin", a_tm_im_scat, g_numRx * g_numTx * g_numBounces * g_numPaths);
  WRITE_BIN("tau_scat.bin", tau_scat, g_numRx * g_numTx * g_numBounces * g_numPaths);
}

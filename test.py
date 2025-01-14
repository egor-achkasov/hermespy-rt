import numpy as np
from rt import compute_paths

# Define inputs
mesh_filepath = __file__[:__file__.rfind('/') + 1] + 'scenes/box.ply'
rx_positions = np.array([[0., 0., 2.5]], dtype=np.float64)
tx_positions = np.array([[0., 0., 2.5]], dtype=np.float64)
rx_velocities = np.array([[0., 0., 0.]], dtype=np.float64)
tx_velocities = np.array([[0., 0., 0.]], dtype=np.float64)
carrier_frequency = 3.0
sampling_rate = 3*1e9
num_rx = 1
num_tx = 1
num_paths = 10000
num_bounces = 3
num_samples = 150

# Call compute_paths
a_te_re, a_te_im, a_tm_re, a_tm_im, tau = compute_paths(
    mesh_filepath,
    rx_positions,
    tx_positions,
    rx_velocities,
    tx_velocities,
    carrier_frequency,
    sampling_rate,
    num_rx,
    num_tx,
    num_paths,
    num_bounces,
    num_samples
)

print(f"TE Gains: {a_te_re + 1.j*a_te_im}")
print(f"TM Gains: {a_tm_re + 1.j*a_tm_im}")
print(f"Delays: {tau}")
print(f"a_te_re min, max: {np.min(a_te_re)}, {np.max(a_te_re)}")
print(f"a_te_im min, max: {np.min(a_te_im)}, {np.max(a_te_im)}")
print(f"a_tm_re min, max: {np.min(a_tm_re)}, {np.max(a_tm_re)}")
print(f"a_tm_im min, max: {np.min(a_tm_im)}, {np.max(a_tm_im)}")

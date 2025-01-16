import numpy as np
import matplotlib.pyplot as plt
from rt import compute_paths

# Define inputs
mesh_filepath = __file__[:__file__.rfind('/') + 1] + 'scenes/box.ply'
rx_positions = np.array([[0., 0., 2.5]], dtype=np.float64)
tx_positions = np.array([[0., 0., 2.5]], dtype=np.float64)
rx_velocities = np.array([[0., 0., 0.]], dtype=np.float64)
tx_velocities = np.array([[0., 0., 0.]], dtype=np.float64)
carrier_frequency = 3.0
num_rx = 1
num_tx = 1
num_paths = 10000
num_bounces = 3

# Call compute_paths
a_te_re_los, a_te_im_los, a_tm_re_los, a_tm_im_los, tau_los, a_te_re_scat, a_te_im_scat, a_tm_re_scat, a_tm_im_scat, tau_scat = compute_paths(
    mesh_filepath,
    rx_positions,
    tx_positions,
    rx_velocities,
    tx_velocities,
    carrier_frequency,
    num_rx,
    num_tx,
    num_paths,
    num_bounces
)

print("#####################")
print("LoS:")
print("#####################")
print(f"shape(a_te_re_los): {a_te_re_los.shape}")
print(f"TE Gains: {a_te_re_los + 1.j*a_te_im_los}")
print(f"TM Gains: {a_tm_re_los + 1.j*a_tm_im_los}")
print(f"Delays: {tau_los}")
print(f"a_te_re_los min, max: {np.min(a_te_re_los)}, {np.max(a_te_re_los)}")
print(f"a_te_im_los min, max: {np.min(a_te_im_los)}, {np.max(a_te_im_los)}")
print(f"a_tm_re_los min, max: {np.min(a_tm_re_los)}, {np.max(a_tm_re_los)}")
print(f"a_tm_im_los min, max: {np.min(a_tm_im_los)}, {np.max(a_tm_im_los)}")

print("\n#####################")
print("Scatter:")
print("#####################")
print(f"shape(a_te_re_scat): {a_te_re_scat.shape}")
print(f"TE Gains: {a_te_re_scat + 1.j*a_te_im_scat}")
print(f"TM Gains: {a_tm_re_scat + 1.j*a_tm_im_scat}")
print(f"Delays: {tau_scat}")
print(f"a_te_re_scat min, max: {np.min(a_te_re_scat)}, {np.max(a_te_re_scat)}")
print(f"a_te_im_scat min, max: {np.min(a_te_im_scat)}, {np.max(a_te_im_scat)}")
print(f"a_tm_re_scat min, max: {np.min(a_tm_re_scat)}, {np.max(a_tm_re_scat)}")
print(f"a_tm_im_scat min, max: {np.min(a_tm_im_scat)}, {np.max(a_tm_im_scat)}")

assert a_te_re_los.shape == a_te_im_los.shape == a_tm_re_los.shape == a_tm_im_los.shape == tau_los.shape
assert a_te_re_scat.shape == a_te_im_scat.shape == a_tm_re_scat.shape == a_tm_im_scat.shape == tau_scat.shape

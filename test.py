import numpy as np
import matplotlib.pyplot as plt
from rt import compute_paths

# Define inputs
mesh_filepath = __file__[:__file__.rfind('/') + 1] + 'scenes/box.ply'
rx_positions = np.array([[0., 0., 2.5]], dtype=np.float64)
tx_positions = np.array([[0., 0., 2.51]], dtype=np.float64)
rx_velocities = np.array([[0., 0., 0.]], dtype=np.float64)
tx_velocities = np.array([[0., 0., 0.]], dtype=np.float64)
carrier_frequency = 3.0
num_rx = 1
num_tx = 1
num_paths = 10000
num_bounces = 3

# Call compute_paths
output = compute_paths(
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
directions_los = output[0]
a_te_los = output[1] + 1.j*output[2]
a_tm_los = output[3] + 1.j*output[4]
tau_los = output[5]
directions_scat = output[6]
a_te_scat = output[7] + 1.j*output[8]
a_tm_scat = output[9] + 1.j*output[10]
tau_scat = output[11]

print("#####################")
print("LoS:")
print("#####################")
print(f"shape(directions_los): {directions_los.shape}")
print(f"shape(a_te_los): {a_te_los.shape}")
print(f"Delays: {tau_los}")
print(f"a_te_los.real min, max: {np.min(a_te_los.real)}, {np.max(a_te_los.real)}")
print(f"a_te_los.imag min, max: {np.min(a_te_los.imag)}, {np.max(a_te_los.imag)}")
print(f"a_tm_los.real min, max: {np.min(a_tm_los.real)}, {np.max(a_tm_los.real)}")
print(f"a_tm_los.imag min, max: {np.min(a_tm_los.imag)}, {np.max(a_tm_los.imag)}")

print("\n#####################")
print("Scatter:")
print("#####################")
print(f"shape(directions_scat): {directions_scat.shape}")
print(f"shape(a_te_scat): {a_te_scat.shape}")
print(f"Delays: {tau_scat}")
print(f"a_te_scat.real min, max: {np.min(a_te_scat.real)}, {np.max(a_te_scat.real)}")
print(f"a_te_scat.imag min, max: {np.min(a_te_scat.imag)}, {np.max(a_te_scat.imag)}")
print(f"a_tm_scat.real min, max: {np.min(a_tm_scat.real)}, {np.max(a_tm_scat.real)}")
print(f"a_tm_scat.imag min, max: {np.min(a_tm_scat.imag)}, {np.max(a_tm_scat.imag)}")

# Asssert shapes
assert a_te_los.shape == a_tm_los.shape == tau_los.shape == (num_rx, num_tx)
assert a_te_scat.shape == a_tm_scat.shape == tau_scat.shape == (num_bounces, num_rx, num_tx, num_paths)

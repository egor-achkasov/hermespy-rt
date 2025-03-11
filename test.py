import numpy as np
from rt import compute_paths

# Define inputs
mesh_filepath = __file__[:__file__.rfind('/') + 1] + 'scenes/simple_reflector.hrt'
rx_positions = np.array([[0., 0., .15]], dtype=np.float64)
tx_positions = np.array([[0., 0., .151]], dtype=np.float64)
rx_velocities = np.array([[0., 0., 0.]], dtype=np.float64)
tx_velocities = np.array([[0., 0., 0.]], dtype=np.float64)
carrier_frequency = 3.0
num_rx = 1
num_tx = 1
num_paths = 10000
num_bounces = 3

# Call compute_paths
los, scatter = compute_paths(
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
print(f"shape(los.directions_rx): {los.directions_rx.shape}")
print(f"shape(los.directions_tx): {los.directions_tx.shape}")
print(f"shape(los.a_te): {los.a_te.shape}")
print(f"shape(los.a_tm): {los.a_tm.shape}")
print(f"Delays: {los.tau}")
print(f"los.a_te.real min, max: {np.min(los.a_te.real)}, {np.max(los.a_te.real)}")
print(f"los.a_te.imag min, max: {np.min(los.a_te.imag)}, {np.max(los.a_te.imag)}")
print(f"los.a_tm.real min, max: {np.min(los.a_tm.real)}, {np.max(los.a_tm.real)}")
print(f"los.a_tm.imag min, max: {np.min(los.a_tm.imag)}, {np.max(los.a_tm.imag)}")

print("\n#####################")
print("Scatter:")
print("#####################")
print(f"shape(scatter.directions_rx): {scatter.directions_rx.shape}")
print(f"shape(scatter.directions_tx): {scatter.directions_tx.shape}")
print(f"shape(scatter.a_te): {scatter.a_te.shape}")
print(f"shape(scatter.a_tm): {scatter.a_tm.shape}")
print(f"Delays: {scatter.tau}")
print(f"scatter.a_te.real min, max: {np.min(scatter.a_te.real)}, {np.max(scatter.a_te.real)}")
print(f"scatter.a_te.imag min, max: {np.min(scatter.a_te.imag)}, {np.max(scatter.a_te.imag)}")
print(f"scatter.a_tm.real min, max: {np.min(scatter.a_tm.real)}, {np.max(scatter.a_tm.real)}")
print(f"scatter.a_tm.imag min, max: {np.min(scatter.a_tm.imag)}, {np.max(scatter.a_tm.imag)}")

# Asssert shapes
# LoS
assert los.num_paths == 1
assert (
    (num_rx, num_tx, 1, 3)
    == los.directions_rx.shape
    == los.directions_tx.shape
)
assert (
    (num_rx, num_tx, 1)
    == los.a_te.shape
    == los.a_tm.shape
    == los.tau.shape
    == los.freq_shift.shape
)
# Scatter
assert scatter.num_paths == num_bounces * num_paths
assert (
    (num_rx, num_tx, scatter.num_paths, 3)
    == scatter.directions_rx.shape
    == scatter.directions_tx.shape
)
assert (
    (num_rx, num_tx, scatter.num_paths)
    == scatter.a_te.shape
    == scatter.a_tm.shape
    == scatter.tau.shape
    == scatter.freq_shift.shape
)

import open3d as o3d
import numpy as np
from sys import argv

pcs = []
for arg in argv[1:]:
    f = open(arg, 'rb')
    #a = np.fromfile(f, dtype=np.float32).reshape((-1, 3))
    #a = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(a))
    #pcs.append(a)
    a = np.fromfile(f, dtype=np.float32).reshape((-1, 10000, 3))
    f.close()
    for ai in a:
        ai = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(ai))
        pcs.append(ai)
pcs.append(o3d.io.read_triangle_mesh('scenes/box.ply'))
o3d.visualization.draw(pcs, show_ui=True)

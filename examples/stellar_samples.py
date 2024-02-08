"""Simulate hererogeneus compound by randomly sampling stellars."""
import random

from bubbles.geometry import Rectangle, Stellar
from bubbles.gmsh_utils import mesh, write_geo
from bubbles.topology import Topology

# Define the domain and create a topology
domain = Rectangle(midpoint=(0, 0), width=1, height=1).to_polygon()
topo = Topology(domain)

levels_unique, weights = zip(*[(1, 5), (2, 4), (3, 3), (4, 2)])

# Sample the levels
n_stellars = 200
levels = random.choices(levels_unique, weights=weights, k=n_stellars)
unif_samples = [random.random() for _ in range(3 * n_stellars)]
midpoints = [
    (mx - 0.5, my - 0.5) for mx, my in zip(unif_samples[::3], unif_samples[1::3])
]

min_r = 0.09
radii = [min(min_r, r * 0.5) for r in unif_samples[2::3]]

for mid, rad, lvl in zip(midpoints, radii, levels):
    C = Stellar(midpoint=mid, radius=rad).to_polygon(refs=50)
    topo.add(polygon=C, level=lvl, extend_domain=False)

topo.plot()
write_geo(
    topology=topo,
    file_name="stellar_samples",
    correct_curve_loops=True,
)

import numpy as np
import scipy.linalg as linalg
import simulate_competitive_cuda as sim

C_total = 1.0
t_pts = np.logspace(-2, 4.5, 300)
t_pts = np.insert(t_pts, 0, 0.0)
site_mismatches = [[]]
# Run with single site (M=1) to eliminate competition factor
P_hist, Cf_hist, clv_hist = sim.simulate_gpu(C_total, site_mismatches, t_pts)

comp_on = clv_hist[:, 0]
nat_on = sim.nature_model(C_total, [], t_pts)

max_diff = np.max(np.abs(comp_on - nat_on))
print(f"Max absolute difference: {max_diff}")
print(f"Final difference: {comp_on[-1] - nat_on[-1]}")

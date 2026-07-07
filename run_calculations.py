import os
import numpy as np
import json

def load_parameters():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    eps_path = os.path.join(base_dir, "dashboard_repo", "parameters", "SpCas9_epsilon.txt")
    rates_path = os.path.join(base_dir, "dashboard_repo", "parameters", "SpCas9_forward_rates.txt")
    
    epsilon = np.loadtxt(eps_path)
    forward_rates = np.loadtxt(rates_path)
    return epsilon, forward_rates

def run():
    epsilon, forward_rates = load_parameters()
    
    # 1. On-target energies (all matches, so mismatch_positions = [])
    # In dead_Cas.py:
    # energies = -1 * epsilon[0:21]
    # energies[0] = epsilon[0]
    energies = -1 * epsilon[0:21]
    energies[0] = epsilon[0]
    
    # 2. Cumulative free energy landscape
    F = [0.0] + list(np.cumsum(energies))
    
    # 3. Mismatch penalties
    delta_epsilon = epsilon[21:]
    
    print("=== ON-TARGET FREE ENERGY DIFFERENCES (Delta F_n) ===")
    for n, e in enumerate(energies):
        state_name = "PAM (0)" if n == 0 else f"State {n}"
        print(f"  {state_name:10s}: {e:12.6f} kBT")
        
    print("\n=== CUMULATIVE FREE ENERGIES (F_n) ===")
    for n, f_val in enumerate(F):
        state_name = "Sol (-1)" if n == 0 else ("PAM (0)" if n == 1 else f"State {n-1}")
        print(f"  {state_name:10s}: {f_val:12.6f} kBT")
        
    print("\n=== MISMATCH PENALTIES (delta_epsilon_n) ===")
    for n, de in enumerate(delta_epsilon):
        print(f"  State {n+1:2d}: {de:12.6f} kBT")
        
    print("\n=== FORWARD RATES ===")
    print(f"  k_on (at 1 nM) : {forward_rates[0]:12.6e} s^-1")
    print(f"  k_f            : {forward_rates[1]:12.6f} s^-1")
    print(f"  k_cat          : {forward_rates[-1]:12.6f} s^-1")

    # Output parameters in JSON
    out_data = {
        "F0": epsilon[0],
        "delta_F": energies.tolist(),
        "F": F,
        "delta_epsilon": delta_epsilon.tolist(),
        "k_on_ref": forward_rates[0],
        "k_f": forward_rates[1],
        "k_cat": forward_rates[-1]
    }
    with open("model_values.json", "w") as f:
        json.dump(out_data, f, indent=4)

if __name__ == "__main__":
    run()

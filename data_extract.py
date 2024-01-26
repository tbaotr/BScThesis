import os, math
import numpy as np
from glob import glob
from scipy.spatial import distance

# ========================== DATA ========================== #
dims = [1, 1, 1, 2, 2, 2, 2, 3, 3, 2, 
        2, 2, 2, 3, 3, 5, 5, 10, 10, 20]
fopts = [-200.0, -1.0, -1.0, -200.0, -1.031628453489877, -186.7309088310239, -1.0, -2709.093505572820, -1.0, 2.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
rhos = [0.01, 0.01, 0.01, 0.01, 0.5, 0.5, 0.2, 0.5, 0.2, 0.01,
        0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
ngopts = [2, 5, 1, 4, 2, 18, 36, 81, 216, 12,
          6, 8, 6, 6, 8, 6, 8, 6, 8, 8]
opt_files = ["F1_opt.dat", "F2_opt.dat", "F3_opt.dat", "F4_opt.dat", "F5_opt.dat", "F6_2D_opt.dat", "F7_2D_opt.dat", "F6_3D_opt.dat", "F7_3D_opt.dat", "F8_2D_opt.dat",
             "CF1_M_D2_opt.dat", "CF2_M_D2_opt.dat", "CF3_M_D2_opt.dat", "CF3_M_D3_opt.dat", "CF4_M_D3_opt.dat", "CF3_M_D5_opt.dat", "CF4_M_D5_opt.dat", "CF3_M_D10_opt.dat", "CF4_M_D10_opt.dat", "CF4_M_D20_opt.dat"]
# ========================================================== #
algos = {"0_0":"HVAM", "0_resi":"ReAM", "5_0":"HVES", "5_resi":"ReES"}
prob_ids = [7]
folder_in = "results"
folder_out = "cec2013"
# ========================================================== #

# Input: results/algo/prob/run/_.dat
# Output: cec2013/prob/algo_prob_dims.txt

if not os.path.exists(f"{folder_out}"):
    os.makedirs(f"{folder_out}")

for prob_id in prob_ids:
    if not os.path.exists(f"{folder_out}/p{prob_id}"):
        os.makedirs(f"{folder_out}/p{prob_id}")

    #1. Find FP and the optimal solution based on CEC2013 benchmark data
    gopts = [j.strip().split() for j in open(f"../[02]_data/cec2013_niching_benchmark/data/{opt_files[prob_id-1]}").readlines()]
    c_gopts = np.array([[float(i) for i in opt] for opt in gopts])

    ## Find the minimum difference between all dimensions
    min_value = 1e9
    for i in range(dims[prob_id-1]):   
        
        temp = np.unique(c_gopts[:, i])    # filtered and sorted
        temp = np.concatenate([temp.reshape(-1,1), temp.reshape(-1,1)], 1)
        dm = distance.cdist(temp, temp, lambda u, v: np.absolute((u-v)).sum()/2)
        
        dm = dm.reshape(-1)
        dm = dm[dm != 0]
        if len(dm) > 0:        
            min_value = min(min_value, np.min(dm))
    if min_value == 1e9:
        min_value = 1

    ## Estimate PF
    PF = math.floor(math.log10(min_value)) - 3  # min_value >= 10^n (n is maximum) => n - 2

    ## Find the optimal solution
    encoded_gopts = []
    for gopt in c_gopts:
        t = ""
        for i in gopt:
            t += hex(math.floor(i * pow(10, -PF)))[2:]
        encoded_gopts.append(t)
    best_solution = "_".join(sorted(set(encoded_gopts)))
    
    #2. Algorithms
    for algo in list(algos.keys()):
        file = open(f"{folder_out}/p{prob_id}/{algos[algo]}_F{prob_id}_{dims[prob_id-1]}D.txt", "w")
        file.writelines("Run\tTime\tPR1\tSR1\tSol1\tPR2\tSR2\tSol2\n")

        ntraject = 0
        runs = glob(f"{folder_in}/{algo}/p{prob_id}/run*")
        for run in runs:
            trajectory = []

            populations = glob(run + "/population*_elites.dat")
            for population in populations:
                e_archive = [j.strip().split() for j in open(population).readlines()]
                c_archive = np.array([[float(i) for i in elite] for elite in e_archive])
                if len(c_archive) == 0:
                    continue
                
                encoded_elites = []
                for elite in c_archive:
                    t = ""
                    for i in elite[:dims[prob_id-1]]:
                        t += hex(math.floor(i * pow(10, -PF)))[2:]
                    encoded_elites.append(t)
                
                dist = distance.cdist(c_archive[:, :dims[prob_id-1]], c_gopts, 'euclidean')               
                
                found_gopts, encoded_population = [], []
                for i in range(len(c_archive)):
                    if abs(c_archive[i, dims[prob_id-1]] - fopts[prob_id-1]) <= 1e-5:
                        found_gopts.append(encoded_gopts[np.argmin(dist[i])])
                        encoded_population.append(encoded_gopts[np.argmin(dist[i])])
                    else:
                        encoded_population.append(encoded_elites[i])
                solution = "_".join(sorted(set(encoded_population)))
                pr = float(len(set(found_gopts))) / float(ngopts[prob_id-1])
                sr = float(len(set(found_gopts))) / float(len(encoded_population))
                
                if len(trajectory) >= 1:
                    if solution != trajectory[-1][2] or sr != trajectory[-1][1] or pr != trajectory[-1][0]:
                        trajectory.append([pr, sr, solution])
                else:
                    trajectory.append([pr, sr, solution])
            
            if len(trajectory) == 0:
                continue
            ntraject += 1
            
            if len(trajectory) == 1:
                file.writelines("{} {1} {:.6f} {:.6f} {} {:.6f} {:.6f} {}\n".format(ntraject, trajectory[0][0], trajectory[0][1], trajectory[0][2], trajectory[0][0], trajectory[0][1], trajectory[0][2]))
            else:
                for i in range(1, len(trajectory)):
                    file.writelines("{} {} {:.6f} {:.6f} {} {:.6f} {:.6f} {}\n".format(ntraject, i, trajectory[i-1][0], trajectory[i-1][1], trajectory[i-1][2], trajectory[i][0], trajectory[i][1], trajectory[i][2]))
        file.close()
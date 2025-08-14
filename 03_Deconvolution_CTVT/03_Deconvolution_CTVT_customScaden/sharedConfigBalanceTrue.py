# shared_config.py

import os
import itertools
import random

BASE_PATH = "/Users/nv4/gitClones/01_Deconvolution/03_Deconvolution_CTVT/033_Results/03_Deconvolution_CTVT_customScaden/0_scadenInput/"
BULK = os.path.join(BASE_PATH, "03_Deconvolution_CTVT_customScaden_bulkCountsGoi.txt")
BULK_PURITY = os.path.join(BASE_PATH, "03_Deconvolution_CTVT_customScaden_bulkCountsGoiPurityAndy.txt")
PURITY = "CTVT"

# SEEDS = [123, 1, 512, 14, 2025]
# BALANCED = [True, False]
# THRESHOLDS = [None, 0, 5, 10]
# RANDOMS = [0, 15, 30, 45]
# VARIANCES = [0.1, 0.2, 0.3, 0.4]
# SAMPLES = [1000, 1500, 2000]
# CELLS = [100, 250, 500]
# LEARN_RATES = [0.0001, 0.001]
# STEPS = [500, 1000, 1500]

SEEDS = [14]
BALANCED = [True]
THRESHOLDS = [None, 0.2]
RANDOMS = [0, 0.3]
VARIANCES = [0]
SAMPLES = [500, 1500] # Training improvements leveled off at 1500 according to paper. 
CELLS = [100, 500]
LEARN_RATES = [0.0001, 0.00001]
STEPS = [1000, 5000]
CELLRANDOM = [.50]

# SEEDS = [123]
# BALANCED = [False]
# THRESHOLDS = [0]
# RANDOMS = [0] 
# VARIANCES = [0.1]
# SAMPLES = [50]
# CELLS = [50]
# LEARN_RATES = [0.01]
# STEPS = [10]

def generate_batches():
    param_grid = list(itertools.product(
        BALANCED, THRESHOLDS, RANDOMS, VARIANCES, SAMPLES, CELLS, LEARN_RATES, STEPS
    ))

    # Filter invalid combos early
    # If balance = True, threshold cannot be set. 
    # If balance = False, threshold must be set, but can include 0 -- none = 0 essentially (testing if simulating with some 0 immune cells sampled helps.)
    param_grid = [
        combo for combo in param_grid
        if (combo[0] is True and combo[1] is None) or
        (combo[0] is False and combo[1] is not None)
    ]

    # Randomize the order of parameter sets
    random.seed(789)

    random.shuffle(param_grid)

    subDirs = [os.path.join(BASE_PATH, d) for d in os.listdir(BASE_PATH)
               if os.path.isdir(os.path.join(BASE_PATH, d))]

    # Use only one subdir for now
    subDir = subDirs[0]

    batches = []
    for params in param_grid:
        bal, thr, ran, var, sam, cel, lr, stp = params
        group = [
            (subDir, bal, thr, ran, var, sam, cel, seed, lr, stp)
            for seed in SEEDS
        ]
        batches.append(group)

    return batches
# shared_config.py

import os
import itertools
import random

BASE_PATH = "/Users/nv4/gitClones/01_Deconvolution/03_Deconvolution_CTVT/033_Results/03_Deconvolution_CTVT_customScaden/snCountsMarkerGenes/0_scadenInput/"
BULK = os.path.join("/Users/nv4/gitClones/01_Deconvolution/03_Deconvolution_CTVT/033_Results/03_Deconvolution_CTVT_customScaden/0_scadenInput/03_Deconvolution_CTVT_customScaden_bulkCountsFinal.txt")
BULK_PURITY = os.path.join("/Users/nv4/gitClones/01_Deconvolution/03_Deconvolution_CTVT/033_Results/03_Deconvolution_CTVT_customScaden/0_scadenInput/03_Deconvolution_CTVT_customScaden_bulkCountsGoiPurityNoBadAndy.txt")
# BULK_PURITY = os.path.join("/Users/nv4/gitClones/01_Deconvolution/03_Deconvolution_CTVT/033_Results/03_Deconvolution_CTVT_customScaden/0_scadenInput/03_Deconvolution_CTVT_customScaden_bulkCountsFinal.txt")
PURITY = "CTVT"
# PURITY = None


SEEDS = [14]
BALANCED = [False]
THRESHOLDS = [0]
RANDOMS = [0] # Mixing between different samples. 
VARIANCES = [0]
SAMPLES = [1500] # Training improvements leveled off at 1500 according to paper. 
CELLS = [500]
LEARN_RATES = [0.00001]
STEPS = [5000]#, 10000]
CELLRANDOM = [0, 0.5] 

def generate_batches():
    param_grid = list(itertools.product(
        BALANCED, THRESHOLDS, RANDOMS, VARIANCES, SAMPLES, CELLS, LEARN_RATES, STEPS, CELLRANDOM
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
    print(subDir)

    batches = []
    for params in param_grid:
        bal, thr, ran, var, sam, cel, lr, stp, cran = params
        group = [
            (subDir, bal, thr, ran, var, sam, cel, seed, lr, stp, cran)
            for seed in SEEDS
        ]
        batches.append(group)

    return batches
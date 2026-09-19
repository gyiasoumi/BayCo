import pandas as pd
import numpy as np



def realseq( 
    # * experimental set up *        
    seed=17,  # random seed (change to get datasets reps for testing)    
    n_background=10000, # number of genes in the background
                        # 10,000 as default as a small RNA-seq dataset
    cluster_sizes=(40, 40, 40), #number of genes in each cluster (number of clusters is implicit)
                                # This default is a bit random 
    timepoints=15, # 10 as default as this is what many datasets have 
    reps=3, #number of biological replicates, 3 defualt as biology standard 

    # * underlying structure of the data * 
    cluster_spread = 0.0001, # standard deviation of the Gaussian cluster in the units of latent parameter space
                            # very small default so that default is co-expressed genes are similar
    cluster_separation = 10, # minimum distance between cluster centres expressed in units of cluster_spread
                             # large default so that default is not co-expressed genes are far apart 
    bg_coexp_fraction = 0, # amount of background genes that are co-expressed
                            # 0 as default 

    # * issues caused by sampling *
    noise= 0.01, # this would a noise of 1% of the mean
                # very low noise as default 
    timepoint_jitter_fraction = 0,  # sampling accuracy (how close to the true timepoint were they actually sampled)
                                    # this value is the sd of the gaussian, this means that 68.2% of the samples will be within this range 
                                    # i.e. if = 0.1, then 68.2% of sample would fall within 10% of the actual timpoint 
                                    # where 10% = 10% of the interval gap betwenn timepoints 
                                    # 0 as default 
    coexpression_mode="positive"   # "positive" or "mixed"
                                    # positive as default for simplicity 
):
    
    # keep here so can easily change
    freq_range=(0.5, 3)
    amp_range=(2.0, 10)
    phase_range=(0.25*np.pi, 1.75*np.pi)
    base_range=(20, 40)

    rng = np.random.default_rng(seed)

    def get_sign():
        if coexpression_mode == "mixed":
            return rng.choice([-1, 1])
        return 1

    t = np.linspace(0, 2*np.pi, timepoints) # spacing the timepoints across the timeseries
    dt = t[1] - t[0]   # spacing between consecutive timepoints (needed for the timepoint jitter)

    ## ** MAKING GENES ** ##

    # helper for getting the parameters for each individual gene profile
    def get_real_params(f_unit, p_unit, a_unit):
        f = freq_range[0] + f_unit * (freq_range[1] - freq_range[0])
        p = phase_range[0] + p_unit * (phase_range[1] - phase_range[0])
        a = amp_range[0] + a_unit * (amp_range[1] - amp_range[0])
        return f, p, a

    ## ** ASSIGNING LATENT CLUSTER SPACES ** ##

    # for the network genes, we must make sure the whole cluster spheres do not overlap
    centers = []
    min_dist = cluster_spread * cluster_separation
    margin = 2 * cluster_spread #to keep the centres away from the edge so that the cluster can fit without being clipped by the edge

    attempts = 0
    while len(centers) < len(cluster_sizes):
        candidate = rng.uniform(margin, 1 - margin, size=3) # (f, p, a) in unit space

        if all(np.linalg.norm(candidate - c) >= min_dist for c in centers):
            centers.append(candidate)

        attempts += 1
        if attempts > 100000:
            raise RuntimeError("Could not place clusters with given spread/separation.")
        
    # getting the co-expressed background genes latent cluster centre

    bg_center = rng.uniform(margin, 1 - margin, size=3) # size = 3 for (f_unit, p_unit, a_unit)
    f_bg_center, p_bg_center, a_bg_center = bg_center    

    ## ** GENERATING THE LATENT PROFILES OF EACH GENE ** ##
    latent_params_all = []

    ## ** BACKGROUND GENES ** ##

    # Number of background genes in the co-expressed cluster
    n_bg_cluster = int(n_background * bg_coexp_fraction)
    n_bg_free = n_background - n_bg_cluster

    for i in range(n_bg_cluster):

        offset = rng.normal(0, cluster_spread, size=3) # how far from the centre of a cluster each gene can spread
        sign = get_sign() # if there is pos and neg co-expression 

        f_unit = f_bg_center # no offset - all genes in a cluster have the same frequency 
        p_unit = p_bg_center + offset[1]
        a_unit = a_bg_center + offset[2]

        latent_params_all.append({
            "Gene": f"Gene_{i:05d}",
            "f_unit": f_unit,
            "p_unit": p_unit,
            "a_unit": a_unit,
            "Cluster": "bg_cluster",
            "Coexpressed": True,
            "base": rng.uniform(*base_range),
            "sign": sign
            })

    for i in range(n_bg_free):

        f_unit, p_unit, a_unit = rng.random(3)
        sign = get_sign()

        latent_params_all.append({
            "Gene": f"Gene_{n_bg_cluster + i:05d}",
            "f_unit": f_unit,
            "p_unit": p_unit,
            "a_unit": a_unit,
            "Cluster": None,
            "Coexpressed": False,
            "base": rng.uniform(*base_range),
            "sign": sign
        })

    ## ** CLUSTER GENES ** ##    
    cluster_ids = []
    cluster_metadata = []

    for c_idx, size in enumerate(cluster_sizes):

        f_center, p_center, a_center = centers[c_idx]

        cluster_metadata.append({
            "cluster": c_idx,
            "f_center": f_center,
            "p_center": p_center,
            "a_center": a_center,
            "spread": cluster_spread
        })

        for _ in range(size):

            offset = rng.normal(0, cluster_spread, size=3)
            sign = get_sign()

            f_unit = f_center
            p_unit = p_center + offset[1]
            a_unit = a_center + offset[2]
        
            cluster_ids.append(c_idx)

            gene_index = n_background + len(cluster_ids) - 1

            latent_params_all.append({
                "Gene": f"Gene_{gene_index:05d}",
                "f_unit": f_unit,
                "p_unit": p_unit,
                "a_unit": a_unit,
                "Cluster": c_idx,
                "Coexpressed": True,
                "base": rng.uniform(*base_range),
                "sign": sign
            })

    ## ** SAMPLING OBSERVED EXPRESSION WITH MEASUREMENT NOISE - THE SIMULATED DATASET (expr_df) ** ## 
    n_total = len(latent_params_all)
    gene_names = [f"Gene_{i:05d}" for i in range(n_total)]

    sample_columns = []
    col_names = []

    for tp in range(timepoints):
        true_t = t[tp]

        for r in range(1, reps+1):
            t_shift = true_t + rng.normal(0, timepoint_jitter_fraction * dt)  # replicate-specific timing offset

            # generate each gene profile on the fly from latent_params_all
            values_list = []

            for d in latent_params_all:
                f, p, a = get_real_params(d['f_unit'], d['p_unit'], d['a_unit'])

                signal = d['sign'] * a * np.sin(f * t + p) + d['base']
                value = np.interp(t_shift, t, signal)

                values_list.append(value)

            values = np.array(values_list)

            # add measurement noise to make the replicates 
            values = values + rng.normal(0, noise * np.abs(values), n_total)

            sample_columns.append(values)
            col_names.append(f"ds_{tp+1}_{r}")

    expr_df = pd.DataFrame(np.column_stack(sample_columns), columns=col_names)
    expr_df.insert(0, "Gene", gene_names)
    
    ## ** RECORDING THE TRUE NETWORK (truth_df) AND WHICH GENES ARE IN EACH CLUSTER (annot) ** ##
    goi_names = gene_names[n_background:] #continue the goi names from the point where the background names finished 
    annot = pd.DataFrame({"Gene": goi_names, "Cluster": cluster_ids})

    c_arr = np.array(cluster_ids)
    truth_df = pd.DataFrame((c_arr[:,None] == c_arr[None,:]).astype(float),
                            index=goi_names, columns=goi_names)
    np.fill_diagonal(truth_df.values, np.nan) # NaNing all self edges

    ## ** saving the input parameters ** ##
    simulation_metadata = {
        "n_background": n_background,
        "cluster_sizes": cluster_sizes,
        "timepoints": timepoints,
        "reps": reps,
        "noise": noise,
        "seed": seed,
        "cluster_spread": cluster_spread,
        "cluster_separation": cluster_separation,
        "bg_coexp_fraction": bg_coexp_fraction,
        "timepoint_jitter_fraction": timepoint_jitter_fraction
    }

    latent_params_df = pd.DataFrame(latent_params_all)

    return expr_df, annot, truth_df, cluster_metadata, simulation_metadata, latent_params_df


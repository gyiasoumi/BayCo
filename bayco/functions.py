# This script includes the functions needed to go from the raw normalised read counts to a network. 

# There are two types of networks which can be built.
## a) networks where the co-expression relationship is compared between every input gene 
## b) networks where the co-exoression relationships are only assessed between pre-specified TF and target genes 


# User Inputs: 

# 1a) a list of strings containing all the genes to be compared 
# 1b) two lists of strings, one containing the TF genes and the other containing the target genes 

# 2) expression data in the following format:
# a single dataset per input file, if multiple datasets are being used, they should be input separately 
## rows = Genes (column can be called "Gene", "gene", "Transcript", or "transcript")
## cols = samples, column names MUST be in the following format X1_X2_X3. where:
### X1 is a dataset identifier string. this should be unique for each dataset, if only using one dataset this still needs to be included
### X2 is a sample identifier string (of sample points e.g. tissue types or treatments), or a number if timepoints. This should be the same across replicates
### X3 is the replicate identifier number
## in this way, every sample has its own unique X1_X2_X3 ID 


# The following functions are for the user to call:

## To build an all-pairs network:
### 1. df_to_prepareddf(), 2. df_to_BFs() (and optionally for testing/troubleshooting df_to_BFs_all_results()), 3. if multiple datasets combine_networks()
## to builf a TF-target pairs network: 
### 1. make_tf_target_pairs() 2. df_to_BFs_tf_vs_targets()

# the following packages need to be installed in an environment 
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.stats import norm, gaussian_kde
from typing import List, Dict, Tuple, Union, Optional
from itertools import product, combinations

#########################################################################
### Preparing the data ###
#########################################################################

################################ HELPERS ################################

# filter for a minimum gene expression
## fiter is the mean across all samples, default is 1 if not specificied by the user 
def filter_for_expression(wide_data: pd.DataFrame, expression_filter: float = 1.0) -> pd.DataFrame:
    wide_data = wide_data.rename(columns={"Transcript": "Gene", "gene": "Gene", "transcript": "Gene"})
    means = wide_data.loc[:, wide_data.columns != 'Gene'].mean(axis=1)
    return wide_data[means > expression_filter].copy()


# converts the wide input dataset into a long format and calculates mean across replicates
def make_long_calc_means(wide_data_filtered: pd.DataFrame) -> pd.DataFrame:

    long_data = pd.melt(wide_data_filtered, id_vars=['Gene'], var_name='Sample', value_name='expr')
    long_data[['ds', 'id', 'rep']] = long_data['Sample'].str.split('_', expand=True)
    
    long_means_sd = long_data.groupby(['Gene', 'id']).agg(
        mean_expr=('expr', 'mean'), 
        sd_expr=('expr', 'std')
    ).reset_index()
    
    long_means_sd['pseudo_sigma'] = np.sqrt(long_means_sd['mean_expr'])
    return long_means_sd

# Z-score normalisation
## Z-score helpers 
def z_mean(series: pd.Series) -> pd.Series:
    return (series - series.mean()) / series.std(ddof=0)

def z_sd(series: pd.Series, mean_std: float) -> pd.Series:
    return series / mean_std

## Z-score normalise the means and sd 
def normalise_means_sd(df: pd.DataFrame) -> pd.DataFrame:

    normalised_df = df.copy()
    normalised_df['z_mean_expr'] = normalised_df.groupby('Gene')['mean_expr'].transform(z_mean)
    
    def normalise_group_sd(g: pd.DataFrame) -> pd.Series:
        return z_sd(g['sd_expr'], g['mean_expr'].std(ddof=0))
        
    normalised_df['z_sd_expr'] = normalised_df.groupby('Gene').apply(
        normalise_group_sd, include_groups=False
    ).reset_index(level=0, drop=True)

    return normalised_df


#########################################################################
# Network Functions of network type a (all pairs) #
#########################################################################

# splitting the prepared dataset into the genes of interest and everything else 
def splitting_df(gene_list: List[str], all_genes_data: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:

    mask = all_genes_data["Gene"].isin(gene_list)
    goi_raw_wide = all_genes_data[mask]
    random_raw_wide = all_genes_data[~mask]
    return goi_raw_wide, random_raw_wide

# calculating the distances for the genes of interest pairs 
def calc_goi_distance_matrix(normalized_df_zscore: pd.DataFrame, mode: str = 'positive') -> pd.DataFrame:

    pivot_df = normalized_df_zscore.pivot(index='id', columns='Gene', values='z_mean_expr')
    genes = pivot_df.columns
    data_matrix = pivot_df.values.T  # (Genes, Timepoints)

    # just the actual distance (e.g. assessing for just positive co-expression)
    diff = data_matrix[:, np.newaxis, :] - data_matrix[np.newaxis, :, :]
    dist_pos = np.sqrt(np.sum(diff**2, axis=2))

    if mode == 'mixed':
        # actual or inverse distance(e.g. assessing for positive and negative co-expression)
        sum_diff = data_matrix[:, np.newaxis, :] + data_matrix[np.newaxis, :, :]
        dist_neg = np.sqrt(np.sum(sum_diff**2, axis=2))
        final_dist = np.minimum(dist_pos, dist_neg)
    elif mode == 'negative':
        # just inverse distance (e.g. assessing for just negative co-expression)
        sum_diff = data_matrix[:, np.newaxis, :] + data_matrix[np.newaxis, :, :]
        final_dist = np.sqrt(np.sum(sum_diff**2, axis=2))
    else:
        final_dist = dist_pos

    final_dist += 1e-10 # to prevent 0 values 
    np.fill_diagonal(final_dist, np.nan)
    return pd.DataFrame(final_dist, index=genes, columns=genes)


# calculating the distances between random background gene pairs ready for the H0 distribution 
def calc_random_distances(random_raw_wide: pd.DataFrame, number_of_pairs: int, mode: str = 'positive') -> pd.DataFrame:

    # pivot_df structure: Index=Timepoints, Columns=Genes
    pivot_df = random_raw_wide.pivot(index='id', columns='Gene', values='z_mean_expr')
    data_matrix = pivot_df.values  # Shape: (n_timepoints, n_genes)
    num_genes = data_matrix.shape[1]
    
    # randomly select indices for all the gene pairs
    # generate slightly more pairs than needed to account for self-matches (p1==p2)
    buffer_size = int(number_of_pairs * 1.1)
    idx1 = np.random.randint(0, num_genes, size=buffer_size)
    idx2 = np.random.randint(0, num_genes, size=buffer_size)
    
    # Filter out self-matches (idx1 == idx2)
    valid_mask = idx1 != idx2
    idx1 = idx1[valid_mask][:number_of_pairs]
    idx2 = idx2[valid_mask][:number_of_pairs]
    
    # Extract columns
    g1_profiles = data_matrix[:, idx1]
    g2_profiles = data_matrix[:, idx2]
    
    # Calculate distances based on Mode
    # actual (positive) distance 
    diff = g1_profiles - g2_profiles
    dist_pos = np.sqrt(np.sum(diff**2, axis=0)) # Sum down the timepoint axis
    
    if mode == 'mixed':
        # actual and inverse distance (mixed)
        sums = g1_profiles + g2_profiles
        dist_neg = np.sqrt(np.sum(sums**2, axis=0))
        # Element-wise minimum
        final_dist = np.minimum(dist_pos, dist_neg)
        
    elif mode == 'negative':
        # inverse distance (negative)
        sums = g1_profiles + g2_profiles
        final_dist = np.sqrt(np.sum(sums**2, axis=0))
        
    else: # 'positive'
        final_dist = dist_pos

    # 5. Format Output
    # Create DataFrame with gene names and distances
    gene_names = pivot_df.columns
    results = pd.DataFrame({
        'Gene1': gene_names[idx1],
        'Gene2': gene_names[idx2],
        'distance': final_dist + 1e-10 # Adding the little bit to prevent 0 values 
    })
    
    return results

# mini wrapper 
def get_distance_dfs(goi_raw_wide: pd.DataFrame, random_raw_wide: pd.DataFrame, 
                     number_of_pairs: int, mode: str = 'positive') -> Tuple[pd.DataFrame, pd.DataFrame]:

    goi_dists = calc_goi_distance_matrix(goi_raw_wide, mode=mode)
    
    rand_dists = calc_random_distances(random_raw_wide, number_of_pairs, mode=mode)
    
    return goi_dists, rand_dists

#########################################################################
# H0 Calculation #
#########################################################################

def fit_gaussian_kde(rand_pair_dist_df: pd.DataFrame) -> gaussian_kde:
    
    data = rand_pair_dist_df['distance'].dropna().values
    return gaussian_kde(data)

def all_get_H0_y_from_gaussian_kde(rand_pair_dist_df: pd.DataFrame, goi_dist_df: pd.DataFrame) -> pd.DataFrame:

    kde = fit_gaussian_kde(rand_pair_dist_df)
    arr = goi_dist_df.values
    y = np.full_like(arr, np.nan, dtype=float)
    mask = ~np.isnan(arr)
    y[mask] = kde(arr[mask])
    return pd.DataFrame(y, index=goi_dist_df.index, columns=goi_dist_df.columns)


#########################################################################
# H1 Calculation #
#########################################################################


def H1_sigma(goi_raw_wide: pd.DataFrame) -> pd.DataFrame:

    pivot_df_sd = goi_raw_wide.pivot(index='id', columns='Gene', values='z_sd_expr')
    ss_sd_vec = np.sum(pivot_df_sd.values**2, axis=0)
    sigma_matrix = np.sqrt(ss_sd_vec[:, np.newaxis] + ss_sd_vec[np.newaxis, :]) + 1e-10
    np.fill_diagonal(sigma_matrix, np.nan)
    return pd.DataFrame(sigma_matrix, index=pivot_df_sd.columns, columns=pivot_df_sd.columns)

def all_H1_from_gaus(goi_raw_wide: pd.DataFrame, goi_dist_df: pd.DataFrame) -> pd.DataFrame:

    H1_sigma_values = H1_sigma(goi_raw_wide)
    mask = (goi_dist_df >= 0).values
    H1_y_array = np.zeros_like(goi_dist_df.values, dtype=float)
    H1_y_array[mask] = 2 * norm.pdf(goi_dist_df.values[mask], loc=0, scale=H1_sigma_values.values[mask])
    H1_y_array[np.isnan(goi_dist_df.values)] = np.nan
    return pd.DataFrame(H1_y_array, index=goi_dist_df.index, columns=goi_dist_df.columns)


#########################################################################
# Bayes Factor Calculation #
#########################################################################

def calc_bayes_matrix(y_H0_df: pd.DataFrame, y_H1_df: pd.DataFrame) -> pd.DataFrame:

    return np.log10(y_H1_df / y_H0_df)


#########################################################################
# Execution Wrappers for network type a (all gene-gene pairs) #
#########################################################################

################## WRAPPER FOR USER to prepare the data #################

def df_to_prepareddf(starting_data: pd.DataFrame, expression_filter: float) -> pd.DataFrame:

    filtered_data = filter_for_expression(starting_data, expression_filter)
    long_genes = make_long_calc_means(filtered_data)

    return normalise_means_sd(long_genes)


################## WRAPPER FOR USER make the network ####################
def df_to_BFs_allgenes(gene_list: List[str], all_genes_data: pd.DataFrame, number_of_pairs: int, mode: str = 'positive') -> pd.DataFrame:
    
    goi_raw_wide, random_raw_wide = splitting_df(gene_list, all_genes_data)
    goi_dist_df, rand_pair_dist_df = get_distance_dfs(goi_raw_wide, random_raw_wide, number_of_pairs, mode=mode)
    y_H0_df = all_get_H0_y_from_gaussian_kde(rand_pair_dist_df, goi_dist_df)
    y_H1_df = all_H1_from_gaus(goi_raw_wide, goi_dist_df)
    return calc_bayes_matrix(y_H0_df, y_H1_df)

# OR to get all the outputs 

def df_to_BFs_allgenes_all_results(gene_list: List[str], all_genes_data: pd.DataFrame, number_of_pairs: int, mode: str = 'positive') -> Dict[str, pd.DataFrame]: 
   
    goi_raw_wide, random_raw_wide = splitting_df(gene_list, all_genes_data)
    goi_dist_df, rand_pair_dist_df = get_distance_dfs(goi_raw_wide, random_raw_wide, number_of_pairs, mode=mode)
    y_H0_df = all_get_H0_y_from_gaussian_kde(rand_pair_dist_df, goi_dist_df)
    y_H1_df = all_H1_from_gaus(goi_raw_wide, goi_dist_df)
    return { 
        'goi_dist_df': goi_dist_df,
        'H0_y_df': y_H0_df, 
        'H1_y_df': y_H1_df, 
        'log10_BF_df': calc_bayes_matrix(y_H0_df, y_H1_df) 
    }

#########################################################################
# Combining Networks #
#########################################################################

# takes in a list of networks built from each df 
# combines networks, if a gene pair is missing in a dataset then it gets a value of 0 for that dataset
# just one function so the user uses this directly

def combine_networks_allgenes(networks: List[pd.DataFrame]) -> pd.DataFrame:

    if not networks:
        raise ValueError("Input list of networks is empty.")

    combined_df = networks[0].copy()

    for df in networks[1:]:
        
        combined_df = combined_df.add(df, fill_value=0)

    np.fill_diagonal(combined_df.values, np.nan)
    return combined_df


#########################################################################
# Instead doing network type b (TF-Target pairs) #
#########################################################################

# make a dictionary of the TF-Target pairs 
def make_tf_target_pairs(tf_list: List[str], target_list: List[str]) -> List[Tuple[str, str]]:
    return [(tf, target) for tf, target in product(tf_list, target_list) if tf != target]

# preparing the data and making the network 
def df_to_BFs_tf_vs_targets(tf_list: List[str], target_list: List[str], all_genes_data: pd.DataFrame, number_of_pairs: int, mode: str = "positive") -> pd.DataFrame:

    # subset to only genes we need
    gene_list = list(set(tf_list) | set(target_list))

    goi_raw_wide, random_raw_wide = splitting_df(gene_list, all_genes_data)

    # calculate random distances (H0 background)
    rand_pair_dist_df = calc_random_distances(random_raw_wide, number_of_pairs, mode=mode)
    kde = fit_gaussian_kde(rand_pair_dist_df)

    # sigma matrix for H1
    sigma_df = H1_sigma(goi_raw_wide)

    # pivot for mean profiles of the TF and target genes
    pivot_df = goi_raw_wide.pivot(index="id", columns="Gene", values="z_mean_expr")

    results = []
    pairs = make_tf_target_pairs(tf_list, target_list)

    for tf, target in pairs:

        if tf not in pivot_df.columns or target not in pivot_df.columns: # need this becasue they could have been filtered out earlier e.g. if are not expressed
            continue

        v1 = pivot_df[tf].values
        v2 = pivot_df[target].values

        # calculating the distances between gene pairs 
        dist_pos = np.sqrt(np.sum((v1 - v2) ** 2))

        if mode == "mixed":
            dist_neg = np.sqrt(np.sum((v1 + v2) ** 2))
            dist = min(dist_pos, dist_neg)
        elif mode == "negative":
            dist = np.sqrt(np.sum((v1 + v2) ** 2))
        else:
            dist = dist_pos

        dist = dist + 1e-10

        # H0
        H0_y = kde(dist)[0]

        # H1
        if tf in sigma_df.index and target in sigma_df.columns:
            sigma = sigma_df.loc[tf, target]
        else:
            sigma = np.nan

        H1_y = 2 * norm.pdf(dist, loc=0, scale=sigma)

        log10_BF = np.log10(H1_y / H0_y)

        results.append((tf, target, dist, log10_BF))
    return pd.DataFrame(results, columns=["TF", "Target", "distance", "log10_BF"])


def combine_tf_target_networks(networks: List[pd.DataFrame]) -> pd.DataFrame:

    if not networks:
        raise ValueError("Input list of networks is empty.")

    combined_df = (
        pd.concat(networks, ignore_index=True)
        .groupby(["TF", "Target"], as_index=False)["log10_BF"]
        .sum()
    )

    return combined_df


#########################################################################
# Other miscellaneous functions  #
#########################################################################

# to inspect to plot the KDE 
def plot_kde(rand_pair_dist_df):
    kde = fit_gaussian_kde(rand_pair_dist_df)
    x = np.linspace(rand_pair_dist_df["distance"].min(),
                    rand_pair_dist_df["distance"].max(), 1000)
    plt.plot(x, kde(x))
    plt.show()
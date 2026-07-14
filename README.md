# BayCo

BayCo description...

## Installation

BayCo requires Python 3.10 or higher.

We recommend installing BayCo in a dedicated conda environment:

First, create the environment:
```bash
conda create -n bayco-env python=3.10
conda activate bayco-env
```
Then, clone the repository and install BayCo and its dependencies:
```
git clone https://github.com/gyiasoumi/BayCo.git
cd BayCo
pip install -e .
```

BayCo can then be imported directly into Python or Jupyter notebooks using
```
import bayco
```

## Input Data Requirements

BayCo requires:

1. A list of genes to analyse, or two separate lists of transcription factors and target genes for TF-target networks.
2. A gene expression dataset supplied as a pandas DataFrame, containing normalised read count metrics such as TPM, FPKM

The expression data must be provided in wide format:

- Rows = genes
- Columns = samples
- Values = normalised expression values

The gene identifier column can be named `Gene`, `gene`, `Transcript`, `transcript`

IMPORTANTLY, sample names must follow the format:

```
DatasetID_SampleID_ReplicateID
```
where:

- `DatasetID` identifies the dataset and must be unique between datasets.
- `SampleID` identifies the biological condition, e.g. tissue, treatment, or timepoint.
- `ReplicateID` identifies biological replicates.

Example:

```
Leaf_0_1
Leaf_0_2
Leaf_24_1
Leaf_24_2
...
```

Multiple datasets should be analysed separately and combined after network inference using `combine_networks()`.


## Building Networks

BayCo can generate two types of networks:

1. **All pairs network**

   Calculates relationships between all input genes.

2. **TF-target network**

   Calculates relationships only between predefined transcription factor and target gene pairs.


## For both types of network, the data first needs to be prepared

```python
prepared_data = df_to_prepareddf(starting_data,
                                 expression_filter=1)
```

This filters genes for a minimum expression level, calculates replicate statistics, and performs the required Z-score normalisation.

### 1. Generate an all pairs network

```python
network = df_to_BFs_allgenes(gene_list,
                    prepared_data,
                    number_of_pairs=30000,
                    mode="positive")
```

Arguments:

- `gene_list`: list of genes to include in the network.
- `prepared_data`: output from `df_to_prepareddf()`.
- `number_of_pairs`: number of random gene pairs used to estimate the background distribution.
- `mode`: type of relationship assessed:
  - `"positive"`: positive co-expression only.
  - `"negative"`: negative co-expression only.
  - `"mixed"`: positive and negative co-expression.

Returns a pandas DataFrame containing pairwise log10 Bayes Factors.


## To Inspecting Intermediate Outputs

To access intermediate calculations to assess qualities of the data which influence the inferences:

```python
results = df_to_BFs_all_results_allgenes(gene_list,
                                         prepared_data,
                                         number_of_pairs=30000)
```

Returns:

- gene-pair distance matrix
- null model likelihoods (H0)
- alternative model likelihoods (H1)
- final log10 Bayes Factor network

To collect evidence of gene pair co-expression across datasets, networks generated from multiple datasets can be combined using:

```python
combined_genenetwork_allgenes = combine_networks_allgenes([network_1, network_2])
```
Missing gene pairs between datasets are assigned a value of zero before combination.



## 2. TF-Target Networks

For predefined transcription factor-target relationships:

```python
tf_target_network = df_to_BFs_tf_vs_targets(tf_list,
                                            target_list,
                                            prepared_data,
                                            number_of_pairs=30000,
                                            mode="positive")
```

Inputs:

- `tf_list`: list of transcription factor genes.
- `target_list`: list of target genes.
- `prepared_data`: prepared expression dataset.

Output columns:

| Column | Description |
|---|---|
| TF | Transcription factor gene |
| Target | Target gene |
| distance | Expression profile distance |
| log10_BF | log10 Bayes Factor |


To collect evidence of gene pair co-expression across datasets, networks generated from multiple datasets can be combined using:

```python
combined_tf_target_genenetwork = combine_tf_target_networks([network_1, network_2])
```

## Examples

Example notebooks demonstrating BayCo workflows are provided in:

```
examples/example_all_pairs_network/
```

and

```
examples/example_TF_target_network/
```
These examples are run on simulated RNA-seq datasets we generated using our own RNA-seq simulator RealSeq.

## RealSeq RNA-seq simulator
We designed RealSeq to test and benchmark BayCo to other association metrics. RealSeq generates RNA-seq datasets with user guide properties such as the number of sample point, number of biological replicates per sample points, 

- the total number of genes in the dataset
- the number of co-expressed genes
- the number co-expression clusters the co-expressed genes are divided into
- the tightness of these co-expressed clusters
- the distance between co-expressed cluster
- the amount non-specific background co-expression in the dataset
- the number of sample points
- the accuracy of the sampling at a timepoint
- the number of replicates per sample points
- the inter-replicate noise. 

RealSeq is available in :

```
RealSeq_data_simulator/
```

## Citation

If you use BayCo in your research, please cite:

[Publication information will be added following publication]


## License

BayCo is released under the MIT License.
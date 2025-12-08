# align_and_tree
This is a simple reproducible nextflow pipeline that will make alignments (using MAFFT) and phylogenetic trees (using iqtree) from one or more sets of nucleotide or protein sequences.

It creates alignments using MAFFT, using default parametes, and then creates trees from the resulting alignments using iqtree (iqtree parameters are [specified in the main nextflow file](https://github.com/stenglein-lab/align_and_tree/blob/eb3cc3eda87509c9f6fd8094fdf50d792c0ff5d3/align_and_tree.nf#L124)).  

This workflow should work for either nucleotide or protein sequences, because neither MAFFT nor iqtree are run in a way that specifies sequence type.  

### Important caveats

- It is the user's responsibility to manually inspect and validate sequence alignments produced by MAFFT. The workflow assumes they are good and just proceeds with tree inference.  The pipeline is not (currently) designed to make trees starting with alignments, but it could be modified to do this.
- The pipeline uses basic default settings for both MAFFT and iqtree.  The pipeline could be modified by editing the code in [the main nextflow file](main.nf). 

### Running the pipeline

A basic command line to run this pipeline would be (assuming the existence of input file of unaligned sequences named `my_sequences.fasta`):

```
nextflow run stenglein-lab/align_and_tree -profile singularity --fasta my_sequences.fasta
```

See [here](https://github.com/stenglein-lab/general_pipeline_instructions) for more information on how to run this pipeline.

#### Multiple input files

It is possible to input multiple fasta files to the pipeline, each of which will produce a separate alignment and tree.  This can be done using a wildcard in the params.fasta value, for instance:

```
nextflow run stenglein-lab/align_and_tree -profile singularity --fasta "input/*.fasta"
```

Note that the wildcard containing params.fasta value must be enclosed in quotes, as [described here](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) and [here](https://www.nextflow.io/docs/latest/cli.html#pipeline-parameters)


### Running test datasets

The pipeline includes a couple small test datasets: sets of influenza A virus segment 4 (HA) [nucleotide](test/influenza_virus_A_HA_nucleotide.fasta) and [protein](test/influenza_virus_A_HA_protein.fasta) RefSeq sequences: N=7 each.  These can be used to confirm that you have correct installations of the software needed to run the pipeline (namely singularity and nextflow).  To run the test datasets, use the [run_test](./run_test) or [run_test_github](./run_test_github) scripts.  [See here](https://github.com/stenglein-lab/general_pipeline_instructions) for more information on dependencies and running test datasets.  These test datasets take ~5 minutes to run on a typical linux multicore server.

### Output

The main outputs will be placed in `results/alignments` and `results/trees` sub-directories.  These include:

| File(s) | Contents | 
| :--- | :--- |
|`results/alignments/*` |  multiple sequence alignment(s) |
|`results/alignments/*.log` |  iqtree log file(s) |
|`results/alignments/*.treefile` |  iqtree .treefile, containing the maximum likelihood tree with support values from both aLRT and bootstrapping.  The support values are in the form: SH-aLRT/UFBoot, where SH-aLRT is the support value from Shimodaira–Hasegawa approximate likelihood ratio testing and the UFBoot is the support value from ultrafast bootstrapping. |
|`results/alignments/*.contree` | iqtree .contree, containing the boostrapping consensus tree and associated support values. |
|`results/alignments/*.iqtree` | file containing log information, various trees in newick and text formats, etc. |

The output directory name can be [overriden](https://www.nextflow.io/docs/latest/workflow.html#publishing-files) using the `-output-dir` command line parameter.

### Optional sequence collapsing prior to alignment

You can optionally collapse sequences prior to alignment using [cd-hit](https://pubmed.ncbi.nlm.nih.gov/23060610/) or [cd-hit-est](https://github.com/weizhongli/cdhit/wiki).  This produces a subset of representatives sequences that each represent clusters of sequences that share some level of pairwise identity. The representative sequences will be aligned.

To do this, specify `--collapse_sequences` on the nextflow command line.  You will also need to specify whether the input fasta files contain nucleotide or protein sequences so the pipeline knows to use cd-hit or cd-hit-est by specifying one of `--sequence_type nt` or `--sequence_type aa`.  You can also specify a pairwide identity threshold for collapsing (input to the cd-hit -c parameter).  Do this using the `--collapse_threshold` parameter (default = 0.98).  

For example:

```
# align and make tree from sequences collapsed at a 99% pariwse identity level
nextflow run stenglein-lab/align_and_tree \
  -profile singularity \
  --fasta my_sequences.fasta \
  --collapse_sequences \
  --sequence_type nt \
  --collapse_threshold 0.99
```

### Dependencies

This workflow uses - and requires - nextflow and singularity to handle dependencies (i.e. MAFFT and iqtree).  See [here](https://github.com/stenglein-lab/general_pipeline_instructions?tab=readme-ov-file#Software-dependencies) for more information on the associated requirements.

### Citations

If you use this tool you should cite the tools it uses, including:

- [MAFFT](https://doi.org/10.1093/molbev/mst010): see also references [listed here](https://mafft.cbrc.jp/alignment/software/)
- [iqtree](https://iqtree.github.io/doc/Home#how-to-cite-iq-tree)
- [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)
- [Singularity](https://doi.org/10.1371/journal.pone.0177459)
- [cd-hit](https://pubmed.ncbi.nlm.nih.gov/23060610/)



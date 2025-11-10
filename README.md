# align_and_tree
This is a simple reproducible nextflow pipeline that will make alignments (using MAFFT) and phylogenetic trees (using iqtree) from one or more sets of nucleotide or protein sequences.

It creates alignments using MAFFT, using default parametes, and then creates trees from the resulting alignments using iqtree (iqtree parameters are [specified in the main nextflow file](https://github.com/stenglein-lab/align_and_tree/blob/eb3cc3eda87509c9f6fd8094fdf50d792c0ff5d3/align_and_tree.nf#L124)).  

This workflow should work for either nucleotide or protein sequences, because neither MAFFT nor iqtree are run in a way that specifies sequence type.  

### Important caveats

- It is the user's responsibility to manually inspect and validate sequence alignments produced by MAFFT. The workflow assumes they are good and just proceeds with tree inference.  The pipeline is not (currently) designed to make trees starting with alignments, but it could be modified to do this.
- The pipeline uses basic default settings for both MAFFT and iqtree.  The pipeline could be modified by editing the code in [align_and_tree.nf](align_and_tree.nf). 

### Running the pipeline

A basic command line to run this pipeline would be (assuming the existence of input file of unaligned sequences named `my_sequences.fasta`):

```
run stenglein-lab/align_and_tree -profile singularity --fasta my_sequences.fasta
```

See [here]() for more information on running the pipeline.

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

### Dependencies

This workflow uses - and requires - nextflow and singularity to handle dependencies (i.e. MAFFT and iqtree).  See [here]() for more information on the associated requirements.

### Citations

If you use this tool you should cite the tools it uses, including:

- [MAFFT](https://doi.org/10.1093/molbev/mst010): see also references [listed here](https://mafft.cbrc.jp/alignment/software/)
- [iqtree](https://iqtree.github.io/doc/Home#how-to-cite-iq-tree)
- [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)
- [Singularity](https://doi.org/10.1371/journal.pone.0177459)



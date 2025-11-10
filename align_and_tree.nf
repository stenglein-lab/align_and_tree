/*
 A simple nextflow workflow to align a set of sequences (mafft) 
 and make a tree from the alignment (iqtree3) 

 The input to this pipeline is one or more sets of sequences 
 in fasta format, specified by the --fasta parameter.

 Mark Stenglein 11/3/2025
 */

workflow {

  main:
    // make alignments
    fasta_ch = Channel.fromPath(params.fasta)
    alignment_workflow(fasta_ch)

    // make trees
    tree_workflow(alignment_workflow.out.alignment)

  publish:
    alignment     = alignment_workflow.out.alignment
    contree       = tree_workflow.out.contree
    iqtree        = tree_workflow.out.iqtree  
    treefile      = tree_workflow.out.treefile  
    treelog       = tree_workflow.out.treelog   
}

// define where main output files will go
output {
    alignment {
        path 'alignments'
        mode 'link'
    }
    contree {
        path 'trees'
        mode 'link'
    }
    iqtree {
        path 'trees'
        mode 'link'
    }
    treefile {
        path 'trees'
        mode 'link'
    }
    treelog {
        path 'trees'
        mode 'link'
    }
}


workflow alignment_workflow {
  take:
    fasta_sequences

  main:
    align_sequences(fasta_sequences)

  emit:
    alignment = align_sequences.out.fasta_alignment
}

workflow tree_workflow {
  take:
    alignment

  main:
    build_tree(alignment)

  emit:
    contree  = build_tree.out.contree
    treefile = build_tree.out.treefile
    iqtree   = build_tree.out.iqtree
    treelog  = build_tree.out.log

}

// align one set of sequences using MAFFT
process align_sequences {
  tag "$fasta"
  label 'process_medium'

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/mafft:7.525--h031d066_1':
    'quay.io/biocontainers/mafft:7.525--h031d066_1' }"

  input:
    path fasta

  output:
    path "*mafft*", includeInputs: false , emit: fasta_alignment

  script:
  // insert "mafft" in output filename, right before final extension
  def extension = fasta.extension
  def new_name  = fasta.name.replaceAll(/\Q${extension}\E$/, "mafft.${extension}")
  """
  mafft ${fasta} > ${new_name}
  """
}

// make a tree using iqtree v3
process build_tree {
  tag "$alignment"
  label 'process_medium'

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/iqtree:3.0.1--h503566f_0':
    'quay.io/biocontainers/iqtree:3.0.1--h503566f_0' }"

  input:
    path alignment

  output:
    path "*.contree",   emit: contree
    path "*.iqtree",    emit: iqtree
    path "*.treefile",  emit: treefile
    path "*.log",       emit: log

  script:
  """
    iqtree -s ${alignment} \\
      -m MFP \\
      -B 1000 \\
      -alrt 1000 \\
      -T auto 
  """
}



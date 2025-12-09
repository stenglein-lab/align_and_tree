/*
 A simple nextflow workflow to align a set of sequences (mafft) 
 and make a tree from the alignment (iqtree3) 

 The input to this pipeline is one or more sets of sequences 
 in fasta format, specified by the --fasta parameter.

 Mark Stenglein 11/3/2025
 */


workflow {

  main:
    // the main input fasta file
    fasta_ch = Channel.fromPath(params.fasta)

    // optional step to collapse nearly identical sequences using cd-hit
    optional_collapsed_output = Channel.empty()
    if (params.collapse_sequences) {

      collapse_workflow(fasta_ch)
      optional_collapsed_output = collapse_workflow.out.collapsed_sequences
      fasta_post_ch = collapse_workflow.out.collapsed_sequences

    } else {
       // just use uncollapsed input ch
       fasta_post_ch = fasta_ch
    }

    // make alignments
    alignment_workflow(fasta_post_ch)

    // make trees
    tree_workflow(alignment_workflow.out.alignment)
  publish:
    collapsed     = optional_collapsed_output
    alignment     = alignment_workflow.out.alignment
    contree       = tree_workflow.out.contree
    iqtree        = tree_workflow.out.iqtree  
    treefile      = tree_workflow.out.treefile  
    treelog       = tree_workflow.out.treelog   
}

// define where main output files will go
output {
    collapsed {
        path 'collapsed_sequences'
        mode 'link'
    }
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


workflow collapse_workflow {
  take:
    fasta_sequences

  main:

    // setup and check collapsing-related parameters

    // pairwise identity threshold
    collapse_threshold = Channel.value(params.collapse_threshold)
    if ((!params.collapse_threshold || params.collapse_threshold < 0 || params.collapse_threshold > 1)){
      error("Error: must define a numeric collapse_threshold in the range of [0-1]")
    }

    // check that sequence type is defined for sequence collapsing
    def allowed_seq_types = ['nt', 'aa']
    if (!params.sequence_type || !allowed_seq_types.contains(params.sequence_type)) {
      error("Error: must define sequence_type as one of  " + allowed_seq_types)
    }
    collapse_seq_type  = Channel.value(params.sequence_type)

    collapse_sequences(fasta_sequences, collapse_threshold, collapse_seq_type)

  emit:
    collapsed_sequences = collapse_sequences.out.collapsed_sequences
}

workflow alignment_workflow {
  take:
    fasta_sequences

  main:
    
    align_sequences(fasta_sequences)

    aligned_output = align_sequences.out.fasta_alignment

    if (params.trim_alignments) {
       trim_mode = Channel.value(params.clipkit_trim_mode)
       trim_alignment(align_sequences.out.fasta_alignment, trim_mode)
       aligned_output = trim_alignment.out.trimmed_alignment
    }

  emit:
    alignment = aligned_output
}

workflow tree_workflow {
  take:
    alignment

  main:
    count_fasta(alignment)
    build_tree(count_fasta.out.fasta_with_count)

  emit:
    contree  = build_tree.out.contree
    treefile = build_tree.out.treefile
    iqtree   = build_tree.out.iqtree
    treelog  = build_tree.out.log

}

// align one set of sequences using MAFFT
process collapse_sequences {
  tag "$fasta"
  label 'process_medium'

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/cd-hit:4.8.1--h5ca1c30_13':
    'quay.io/biocontainers/cd-hit:4.8.1--h5ca1c30_13' }"

  input:
    path fasta
    val threshold
    val sequence_type

  output:
    path "*.cdhit.*", includeInputs: false , emit: collapsed_sequences

  script:

  // cd-hit or cd-hit-est?
  def tool = sequence_type == "nt" ? "cd-hit-est" : "cd-hit"

  // insert "cdhit" in output filename, right before final extension
  def extension = fasta.extension
  def new_name  = fasta.name.replaceAll(/\Q${extension}\E$/, "cdhit.${threshold}.${extension}")

  // cd-hit expects memory limits in MB
  mem = task.memory.getMega() 
  """
  # run cd-hit or cd-hit-est
  ${tool} -c ${threshold} -i ${fasta} -o ${new_name} -M ${mem} -T ${task.cpus}
  # get rid of .clstr output file so it doesn't populate output
  rm -f ${new_name}.clstr
  """
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
  // TODO: setup/control memory and cpu usage
  """
  mafft ${fasta} > ${new_name}
  """
}

// Count # of sequences in a fasta file
process count_fasta {
  tag "$fasta"
  label 'process_low'

  input:
    path fasta

  output:
    tuple path ("$fasta"), stdout, emit: fasta_with_count

  script:
  """
  grep -c \\> $fasta | tr -d '\n'
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
    tuple path (alignment), val(num_seqs)

  output:
    path "*.contree",   emit: contree, optional: true
    path "*.iqtree",    emit: iqtree
    path "*.treefile",  emit: treefile
    path "*.log",       emit: log

  script:
  mem       = task.memory.getMega() + "M"
  // only bootstrap if at least 4 sequences
  int int_num_seqs = Integer.parseInt(num_seqs)
  bootstrap = int_num_seqs >= 4 ? "-B 1000" : ""

  """
    iqtree -s ${alignment} \\
      -m MFP \\
      $bootstrap \\
      -alrt 1000 \\
      -T ${task.cpus} \\
      -mem ${mem}
  """
}


// trim alignment with clipkit
process trim_alignment {
  tag "$alignment"
  label 'process_medium'

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/clipkit:2.7.0--pyhdfd78af_0':
    'quay.io/biocontainers/clipkit:2.7.0--pyhdfd78af_0' }"

  input:
    path alignment
    val trim_mode

  output:
    path "*clipkit.*",   emit: trimmed_alignment

  script:

  // TODO: allow 
  def extension = alignment.extension
  def new_name  = alignment.name.replaceAll(/\Q${extension}\E$/, "clipkit.${extension}")
  """
    clipkit \\
      $alignment \\
      -m $trim_mode \\
      -o $new_name \\
      -t ${task.cpus} 
  """
}

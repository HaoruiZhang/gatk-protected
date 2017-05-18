# This workflow is used for running germline CNV on a cohort of germline samples
# Notes: 
#
# -Basic sex genotype tab separated table for homo sapiens must be formatted as follows (See SexGenotypeTableReader.java class for full description):
#    SAMPLE_NAME  SEX_GENOTYPE
#    sample_name_1 SEX_XX
#    sample_name_2 SEX_XY
#    sample_name_3  SEX_XY
#    sample_name_4  SEX_XX
# where sex genotype identifiers must match those in tab-separated contig ploidy annotation table that should be formatted as follows:
#    CONTIG  CLASS  SEX_XX  SEX_XY
#    1         AUTOSOMAL       2          2
#    2         AUTOSOMAL       2          2
#    ...       ...             ...        ...
#    X         ALLOSOMAL       2          0
#    Y         ALLOSOMAL       1          1
#
# - Input file (normal_bams_list) must contain file paths to bam and bam index files separated by tabs in the following format:
#    normal_bam_1 bam_idx_1
#    normal_bam_2 bam_idx_2
#
# - The target file (targets) is required for the WES workflow and should be a TSV file with the column headers:
#    contig    start    stop    name
#   These targets will be padded on both sides by the amount specified by PadTargets.padding (default 250).
#
# - If a target file is not provided, then the WGS workflow will be run instead and the specified value of
#   wgs_bin_size (default 10000) will be used.
#
# - Example invocation:
#    java -jar cromwell.jar run gCNV_cohort_calling_workflow.wdl myParameters.json
#   See gCNV_cohort_calling_workflow.json for a template json file to modify with your own parameters (please save
#   your modified version with a different filename and do not commit to the gatk-protected repository).
################


import "gCNV_single_sample_calling_workflow.wdl" as gCNVCalling

workflow gCNVCohortCallingWorkflow {
  # Workflow input files
  File? targets
  File normal_bams_list
  Array[Array[String]]+ normal_bams = read_tsv(normal_bams_list)
  File ref_fasta
  File ref_fasta_dict
  File ref_fasta_fai
  File sex_genotypes
  File contig_annotations
  String gatk_jar

  # Transition prior table files
  File transition_prior_table
  File? transition_matrix_XX_Y
  File? transition_matrix_XY_X
  File? transition_matrix_XY_Y
  File? transition_matrix_XX_X
  File? transition_matrix_autosomal

  # Model directory and parameters
  String model_path
  Int num_latents

  # Output path
  String output_path

  scatter (normal_bam in normal_bams) {
    call gCNVCalling.gCNVSingleSampleWorkflow as CohortCalling {
      input:
        targets = targets,
        normal_bam = normal_bam[0],
        normal_bam_idx = normal_bam[1],
        ref_fasta = ref_fasta,
        ref_fasta_dict = ref_fasta_dict,
        ref_fasta_fai = ref_fasta_fai,
        sex_genotypes = sex_genotypes,
        contig_annotations = contig_annotations,
        gatk_jar = gatk_jar,
        transition_prior_table = transition_prior_table,
        transition_matrix_XX_Y = transition_matrix_XX_Y,
        transition_matrix_XY_X = transition_matrix_XY_X,
        transition_matrix_XY_Y = transition_matrix_XY_Y,
        transition_matrix_XX_X = transition_matrix_XX_X,
        transition_matrix_autosomal = transition_matrix_autosomal,
        model_path = model_path,
        output_path = output_path,
        num_latents = num_latents
    }
  }
  
  output {
    Array[Array[File]] posterior_files = CohortCalling.posteriors
    Array[Array[File]] segment_files = CohortCalling.segments
  }
  
}

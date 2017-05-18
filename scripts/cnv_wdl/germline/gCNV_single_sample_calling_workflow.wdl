# Subworkflow for running GATK germline CNV on a single BAM. Supports both WGS and WES samples.
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
# - The target file (targets) is required for the WES workflow and should be a TSV file with the column headers:
#    contig    start    stop    name
#   These targets will be padded on both sides by the amount specified by PadTargets.padding (default 250).
#
# - If a target file is not provided, then the WGS workflow will be run instead and the specified value of
#   wgs_bin_size (default 10000) will be used.
#
# - Example invocation:
#    java -jar cromwell.jar run gCNV_single_sample_calling_workflow.wdl myParameters.json
#   See gCNV_single_sample_calling_workflow.json for a template json file to modify with your own parameters (please save
#   your modified version with a different filename and do not commit to the gatk-protected repository).
################

import "cnv_common_tasks.wdl" as CNVTasks

workflow gCNVSingleSampleWorkflow {
  # Workflow input files
  File? targets  
  File normal_bam
  File normal_bam_idx
  File ref_fasta
  File ref_fasta_dict
  File ref_fasta_fai
  File sex_genotypes
  File contig_annotations 
  String gatk_jar

  # Transtion prior table files
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

  # If no target file is input, then do WGS workflow
  Boolean is_wgs = !defined(targets) 
  
  if (!is_wgs) {
    call CNVTasks.PadTargets {
      input:
        targets = targets,
        gatk_jar = gatk_jar
    }
  }

  call CNVTasks.CollectCoverage {
    input:
      padded_targets = PadTargets.padded_targets,
      bam = normal_bam,
      bam_idx = normal_bam_idx,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      ref_fasta_dict = ref_fasta_dict,
      gatk_jar = gatk_jar,
      transform = "RAW"
  }

  call CNVTasks.AnnotateTargets {
    input:
      entity_id = CollectCoverage.entity_id,
      targets = CollectCoverage.coverage,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      ref_fasta_dict = ref_fasta_dict,
      gatk_jar = gatk_jar
  }

  call CNVTasks.CorrectGCBias {
    input:
      entity_id = CollectCoverage.entity_id,
      coverage = CollectCoverage.coverage,
      annotated_targets = AnnotateTargets.annotated_targets,
      gatk_jar = gatk_jar
  }

  call GermlineCNVCaller {
    input:
      coverage = CorrectGCBias.corrected_coverage,
      contig_annotations = contig_annotations,
      sex_genotypes = sex_genotypes,
      transition_prior_table = transition_prior_table,
      transition_matrix_XX_Y = transition_matrix_XX_Y,
      transition_matrix_XY_X = transition_matrix_XY_X,
      transition_matrix_XY_Y = transition_matrix_XY_Y,
      transition_matrix_XX_X = transition_matrix_XX_X,
      transition_matrix_autosomal = transition_matrix_autosomal,
      model_path = model_path,
      num_latents = num_latents, 
      output_path = output_path,
      gatk_jar = gatk_jar
  }
 
  output {
    Array[File] posteriors = GermlineCNVCaller.posteriors
    Array[File] segments = GermlineCNVCaller.segments
  }  
}

# Call CNV events on germline samples using panel of normals
task GermlineCNVCaller {
    File coverage
    File contig_annotations
    File sex_genotypes
    File transition_prior_table
    File? transition_matrix_XX_Y
    File? transition_matrix_XY_X
    File? transition_matrix_XY_Y
    File? transition_matrix_XX_X
    File? transition_matrix_autosomal
    String output_path
    String model_path
    Int num_latents
    File gatk_jar
    Int? mem

    command {
        java -Xmx${default=4 mem}g -Ddtype=double -jar ${gatk_jar} GermlineCNVCaller \
            --input ${coverage} \
            --inputModelPath ${model_path} \
            --contigAnnotationsTable ${contig_annotations} \
            --sexGenotypeTable ${sex_genotypes} \
            --copyNumberTransitionPriorTable ${transition_prior_table} \
            --outputPath ${output_path} \
            --numLatents ${default=5 num_latents} \
            --jobType CALL_ONLY \
            --rddCheckpointing false \
            --disableSpark true
    }

    output {
        Array[File] posteriors = glob("./${output_path}/posteriors_final/*")  
        Array[File] segments = glob("./${output_path}/posteriors_final/segments/*")
    }
}

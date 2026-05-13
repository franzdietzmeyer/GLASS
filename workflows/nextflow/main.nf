nextflow.enable.dsl = 2

/*
 * GLASS Nextflow entry — sole workflow driver.
 * Params are provided via -params-file from glass_config_json.py (see run_glass_nextflow.sh).
 *
 * Rosetta steps use errorStrategy 'ignore' and scripts end with || true so sibling jobs continue;
 * merge steps use filesystem globs so partial successes still produce a merged scorefile when possible.
 *
 * DEBUG: GLASS_NEXTFLOW_DEBUG=1 in the environment (see scripts and launcher).
 *
 * Per-task logs under results/.../logs/: scripts/nextflow_copy_task_logs.sh (NF DSL cannot call top-level def helpers).
 *
 * Resource labels (see workflows/nextflow/nextflow.config and conf/slurm.config):
 *   rosetta_initial_relax — FastRelax replicate (container).
 *   rosetta_nogly         — glycan masking, no_glycans or no_glycans_forced (run_one_position.sh).
 *   rosetta_glycans       — glycan masking, glycans or glycans_forced: parallel workers per position (shared scorefile, MPWOD).
 *   nf_light              — prepare_positions, finalize, merges, analysis.
 */

process INITIAL_RELAX_REPLICATE {
    tag { "initial_relax_replicate_${replicate_id}" }
    label 'rosetta_initial_relax'
    cpus 1
    afterScript { "bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'initial_relax_replicate' '${task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')}'" }

    input:
    val(replicate_id)

    output:
    val(replicate_id), emit: replicate_done

    script:
    def out_dir_rel = new File(params.pdb_path as String).parent ?: '.'
    // replicate_id only labels the Nextflow task; each run uses the same Rosetta command and paths.
    """
    cd '${params.launch_dir}'
    export GLASS_CONFIG_INI='${params.config_ini_abs}'
    bash scripts/run_initial_relax_replicate.sh \\
      '${params.input_pdb_path}' \\
      '${params.config_ini_abs}' \\
      '${out_dir_rel}'
    """
}

process INITIAL_RELAX_FINALIZE {
    tag 'initial_relax_finalize'
    label 'nf_light'
    cpus 1
    afterScript { "bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'initial_relax_finalize' '${task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')}'" }

    input:
    val(replicate_ids)

    output:
    path 'out/relaxed.pdb', emit: relaxed_pdb

    script:
    def out_dir_rel = new File(params.pdb_path as String).parent ?: '.'
    // Nextflow 25+ only allows path outputs inside the task work dir — copy the canonical PDB here for staging.
    """
    WORK_DIR="\$PWD"
    cd '${params.launch_dir}'
    export GLASS_CONFIG_INI='${params.config_ini_abs}'
    bash scripts/run_initial_relax_finalize.sh \\
      '${params.config_ini_abs}' \\
      '${out_dir_rel}' \\
      '${params.pdb_path}'
    mkdir -p "\$WORK_DIR/out"
    cp '${params.pdb_path_abs}' "\$WORK_DIR/out/relaxed.pdb"
    """
}

process PREPARE_POSITIONS {
    tag 'prepare_positions'
    label 'nf_light'
    cpus 1
    afterScript { "bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'prepare_positions' '${task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')}'" }

    input:
    val _ready

    output:
    val(true), emit: done

    script:
    def dbg = params.pipeline_debug ? 'True' : 'False'
    """
    cd '${params.launch_dir}'
    export GLASS_CONFIG_INI='${params.config_ini_abs}'
    export GLASS_INPUT_PDB='${params.pdb_path}'
    bash scripts/run_prepare_positions.sh '${params.positions_dir}' ${dbg} '${params.pdb_path}'
    """
}

process GLYCAN_MASKING_NOGLY {
    label 'rosetta_nogly'
    cpus 1
    errorStrategy 'ignore'
    maxRetries 0
    tag { "mask_${position_id}" }
    afterScript { "bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'glycan_masking_nogly' '${task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')}'" }

    input:
    val position_id

    output:
    val(position_id), emit: done

    script:
    """
    cd '${params.launch_dir}'
    export GLASS_CONFIG_INI='${params.config_ini_abs}'
    bash scripts/run_glycan_masking.sh \\
      '${position_id}' \\
      '${params.pdb_path}' \\
      '${params.config_ini_abs}' \\
      '${params.result_dir}/out_by_position/${position_id}' \\
      || true
    """
}

// Parallel chunk workers: same position may run multiple tasks (batch_id). Each uses chunk nstruct,
// shared per-position scorefile + MPWOD + stagger (run_one_position.sh).
process GLYCAN_MASKING_BATCH {
    label 'rosetta_glycans'
    cpus 1
    errorStrategy 'ignore'
    maxRetries 0
    tag { "mask_${position_id}_b${batch_id}" }
    afterScript { "bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'glycan_masking_batch' '${task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')}'" }

    input:
    tuple val(position_id), val(batch_id)

    output:
    tuple val(position_id), val(batch_id), emit: done

    script:
    """
    cd '${params.launch_dir}'
    export GLASS_CONFIG_INI='${params.config_ini_abs}'
    bash scripts/run_glycan_masking.sh \\
      '${position_id}' \\
      '${params.pdb_path}' \\
      '${params.config_ini_abs}' \\
      '${params.result_dir}/out_by_position/${position_id}' \\
      '${batch_id}' \\
      '${params.glycan_batch_size}' \\
      '${params.nstruct}' \\
      || true
    """
}

process GLYCAN_MERGE_POSITIONS {
    tag 'merge_position_batches'
    label 'nf_light'
    cpus 1
    afterScript { "bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'glycan_merge_positions' '${task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')}'" }

    input:
    val _batch_done

    output:
    val(true), emit: done

    script:
    """
    cd '${params.launch_dir}'
    bash scripts/nextflow_merge_glycan_positions.sh \\
      '${params.launch_dir}' \\
      '${params.workflow_pdb_stem ?: params.pdb_name}' \\
      '${params.result_dir}' \\
      '${params.positions_list}'
    """
}

process GLOBAL_MERGE {
    tag 'merge_scores'
    label 'nf_light'
    cpus 1
    afterScript { "bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'global_merge' '${task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')}'" }

    input:
    val _ready

    output:
    val(true), emit: done

    script:
    """
    cd '${params.launch_dir}'
    bash scripts/nextflow_merge_global.sh \\
      '${params.launch_dir}' \\
      '${params.workflow_pdb_stem ?: params.pdb_name}' \\
      '${params.result_dir}' \\
      '${params.merged_score_path}'
    """
}

process ANALYZE {
    tag 'analyze'
    label 'nf_light'
    cpus 1
    afterScript { "bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'analyze' '${task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')}'" }

    input:
    val _ready

    output:
    val(true), emit: done

    script:
    """
    cd '${params.launch_dir}'
    bash scripts/run_analysis.sh \\
      '${params.config_ini_abs}' \\
      '${params.merged_score_path}' \\
      '${params.pdb_path}' \\
      '${params.workflow_pdb_stem ?: params.pdb_name}' \\
      '${params.analysis_marker}' \\
      '${params.out_dir}'
    """
}

process POSITION_RUN_AUDIT {
    tag 'position_run_audit'
    label 'nf_light'
    cpus 1
    afterScript { "bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'position_run_audit' '${task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')}'" }

    input:
    val _ready

    output:
    val(true), emit: done

    script:
    """
    cd '${params.launch_dir}'
    bash scripts/position_run_audit.sh \\
      '${params.launch_dir}' \\
      '${params.result_dir}' \\
      '${params.positions_list}'
    """
}

// -----------------------------------------------------------------------------
// Entry workflow (single block for compatibility across Nextflow 24+)
// -----------------------------------------------------------------------------

workflow {
    main:
    if (params.initial_relax == true) {
        ch_initial_relax_replicates = channel.from(1..params.initial_relax_nstruct)
        INITIAL_RELAX_REPLICATE(ch_initial_relax_replicates)
        INITIAL_RELAX_FINALIZE(INITIAL_RELAX_REPLICATE.out.replicate_done.collect())
        // PREPARE must run only after the relaxed PDB exists (output path staged above).
        ch_after_relax = INITIAL_RELAX_FINALIZE.out.relaxed_pdb.map { true }
    } else {
        ch_after_relax = channel.of(true)
    }
    PREPARE_POSITIONS(ch_after_relax)
    ch_pos = PREPARE_POSITIONS.out.done
        .map { file("${params.launch_dir}/${params.positions_list}") }
        .splitText()
        .map { it.trim() }
        .filter { it }

    if (params.glycan_model == 'glycans' || params.glycan_model == 'glycans_forced') {
        ch_batches = Channel.from(1..params.n_batches)
        ch_jobs = ch_pos.combine(ch_batches)
        GLYCAN_MASKING_BATCH(ch_jobs)
        GLYCAN_MERGE_POSITIONS(GLYCAN_MASKING_BATCH.out.done.collect())
        GLOBAL_MERGE(GLYCAN_MERGE_POSITIONS.out.done)
        POSITION_RUN_AUDIT(GLOBAL_MERGE.out.done)
        ANALYZE(POSITION_RUN_AUDIT.out.done)
    } else {
        GLYCAN_MASKING_NOGLY(ch_pos)
        GLOBAL_MERGE(GLYCAN_MASKING_NOGLY.out.done.collect())
        POSITION_RUN_AUDIT(GLOBAL_MERGE.out.done)
        ANALYZE(POSITION_RUN_AUDIT.out.done)
    }
}

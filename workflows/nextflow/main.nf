nextflow.enable.dsl = 2

/*
 * GLASS Nextflow entry (additive; does not replace Snakemake).
 * Params are provided via -params-file from glass_config_json.py (see run_glass_nextflow.sh).
 *
 * Rosetta steps use errorStrategy 'ignore' and scripts end with || true so sibling jobs continue;
 * merge steps use filesystem globs so partial successes still produce a merged scorefile when possible.
 *
 * DEBUG: GLASS_NEXTFLOW_DEBUG=1 in the environment (see scripts and launcher).
 *
 * Per-task logs under results/.../logs/: scripts/nextflow_copy_task_logs.sh (NF DSL cannot call top-level def helpers).
 */

process INITIAL_RELAX {
    tag 'initial_relax'
    label 'rosetta'
    cpus { params.local_queue_size }
    afterScript = {
        def safe = task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')
        """
        bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'initial_relax' '${safe}'
        """.stripIndent()
    }

    when:
    params.initial_relax == true

    script:
    // Nextflow 25+ only allows path outputs inside the task work dir — copy the canonical PDB here for staging.
    // run_initial_relax.sh still writes to results/.../initial_relax/<pdb_name>.pdb (repo-relative paths after cd).
    """
    WORK_DIR="\$PWD"
    cd '${params.launch_dir}'
    export GLASS_CONFIG_INI='${params.config_ini_abs}'
    bash scripts/run_initial_relax.sh \\
      '${params.input_pdb_path}' \\
      '${params.config_ini_abs}' \\
      '${params.pdb_path}'
    mkdir -p "\$WORK_DIR/out"
    cp '${params.pdb_path_abs}' "\$WORK_DIR/out/relaxed.pdb"
    """
    output:
    path 'out/relaxed.pdb', emit: relaxed_pdb
}

process PREPARE_POSITIONS {
    tag 'prepare_positions'
    cpus 1
    afterScript = {
        def safe = task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')
        """
        bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'prepare_positions' '${safe}'
        """.stripIndent()
    }

    input:
    val _ready

    script:
    def dbg = params.snakemake_debug ? 'True' : 'False'
    """
    cd '${params.launch_dir}'
    export GLASS_CONFIG_INI='${params.config_ini_abs}'
    export GLASS_INPUT_PDB='${params.pdb_path}'
    bash scripts/run_prepare_positions.sh '${params.positions_dir}' ${dbg} '${params.pdb_path}'
    """
    output:
    val(true), emit: done
}

process GLYCAN_MASKING_NOGLY {
    label 'rosetta'
    cpus 1
    errorStrategy 'ignore'
    maxRetries 0
    afterScript = {
        def safe = task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')
        """
        bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'glycan_masking_nogly' '${safe}'
        """.stripIndent()
    }

    input:
    val position_id

    tag "mask_${position_id}"

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
    output:
    val(position_id), emit: done
}

process GLYCAN_MASKING_BATCH {
    label 'rosetta'
    cpus 1
    errorStrategy 'ignore'
    maxRetries 0
    afterScript = {
        def safe = task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')
        """
        bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'glycan_masking_batch' '${safe}'
        """.stripIndent()
    }

    input:
    tuple val(position_id), val(batch_id)

    tag "mask_${position_id}_b${batch_id}"

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
    output:
    tuple val(position_id), val(batch_id), emit: done
}

process GLYCAN_MERGE_POSITIONS {
    tag 'merge_position_batches'
    cpus 1
    afterScript = {
        def safe = task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')
        """
        bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'glycan_merge_positions' '${safe}'
        """.stripIndent()
    }

    input:
    val _batch_done

    script:
    """
    cd '${params.launch_dir}'
    bash scripts/nextflow_merge_glycan_positions.sh \\
      '${params.launch_dir}' \\
      '${params.pdb_name}' \\
      '${params.result_dir}' \\
      '${params.positions_list}'
    """
    output:
    val(true), emit: done
}

process GLOBAL_MERGE {
    tag 'merge_scores'
    cpus 1
    afterScript = {
        def safe = task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')
        """
        bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'global_merge' '${safe}'
        """.stripIndent()
    }

    input:
    val _ready

    script:
    """
    cd '${params.launch_dir}'
    bash scripts/nextflow_merge_global.sh \\
      '${params.launch_dir}' \\
      '${params.pdb_name}' \\
      '${params.result_dir}' \\
      '${params.merged_score_path}'
    """
    output:
    val(true), emit: done
}

process ANALYZE {
    tag 'analyze'
    cpus 1
    afterScript = {
        def safe = task.name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')
        """
        bash '${params.launch_dir}/scripts/nextflow_copy_task_logs.sh' '${params.launch_dir}' '${params.result_dir}' 'analyze' '${safe}'
        """.stripIndent()
    }

    input:
    val _ready

    script:
    """
    cd '${params.launch_dir}'
    bash scripts/run_analysis.sh \\
      '${params.config_ini_abs}' \\
      '${params.merged_score_path}' \\
      '${params.pdb_path}' \\
      '${params.pdb_name}' \\
      '${params.analysis_marker}' \\
      '${params.out_dir}'
    """
    output:
    val(true), emit: done
}

// -----------------------------------------------------------------------------
// Entry workflow (single block for compatibility across Nextflow 24+)
// -----------------------------------------------------------------------------

workflow {
    main:
    if (params.initial_relax == true) {
        INITIAL_RELAX()
        // PREPARE must run only after the relaxed PDB exists (output path staged above).
        ch_after_relax = INITIAL_RELAX.out.relaxed_pdb.map { true }
    } else {
        ch_after_relax = channel.of(true)
    }
    PREPARE_POSITIONS(ch_after_relax)
    ch_pos = PREPARE_POSITIONS.out.done
        .map { file("${params.launch_dir}/${params.positions_list}") }
        .splitText()
        .map { it.trim() }
        .filter { it }

    if (params.glycan_model == 'glycans') {
        ch_batches = Channel.from(1..params.n_batches)
        ch_jobs = ch_pos.combine(ch_batches)
        GLYCAN_MASKING_BATCH(ch_jobs)
        GLYCAN_MERGE_POSITIONS(GLYCAN_MASKING_BATCH.out.done.collect())
        GLOBAL_MERGE(GLYCAN_MERGE_POSITIONS.out.done)
        ANALYZE(GLOBAL_MERGE.out.done)
    } else {
        GLYCAN_MASKING_NOGLY(ch_pos)
        GLOBAL_MERGE(GLYCAN_MASKING_NOGLY.out.done.collect())
        ANALYZE(GLOBAL_MERGE.out.done)
    }
}

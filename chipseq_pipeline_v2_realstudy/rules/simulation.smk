rule build_prototype_run_table:
    output:
        csv=config["params_table"],
        meta="metadata/prototype_run_table.summary.json"
    params:
        study_ids=lambda wc: " ".join(config["study_ids"]),
        study_depths_csv=lambda wc: config["study_depths_csv"],
        coverage_treat=lambda wc: " ".join(str(value) for value in config["coverage_treat"]),
        coverage_ctrl=lambda wc: " ".join(str(value) for value in config["coverage_ctrl"]),
        seed=lambda wc: " ".join(str(value) for value in config["seed"]),
        fragment_length=lambda wc: " ".join(str(value) for value in config["fragment_length"]),
        read_length=lambda wc: " ".join(str(value) for value in config["read_length"]),
        aligners=lambda wc: " ".join(config["aligners"]),
        peakcallers=lambda wc: " ".join(config["peakcaller_list"]),
        macs2_mode=lambda wc: " ".join(config["macs2_mode"])
    shell:
        """
        python -m scripts.sample_reads_from_intensity \
          --study-ids {params.study_ids} \
          --study-depths-csv {params.study_depths_csv} \
          --observed-treatment-depth 10 \
          --observed-control-depth 4 \
          --coverage-treat {params.coverage_treat} \
          --coverage-ctrl {params.coverage_ctrl} \
          --seed {params.seed} \
          --fragment-length {params.fragment_length} \
          --read-length {params.read_length} \
          --aligners {params.aligners} \
          --peakcaller-list {params.peakcallers} \
          --macs2-mode {params.macs2_mode} \
          --output-csv {output.csv} \
          --output-json {output.meta}
        """

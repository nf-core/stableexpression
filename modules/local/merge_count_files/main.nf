process MERGE_COUNT_FILES {

    // label 'error_retry'
    debug true

    conda "${moduleDir}/environment.yml"

    input:
    val csv_files_str


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    from pathlib import Path

    # changing input value to list of Paths
    csv_files = [Path(file) for file in "$csv_files_str".strip('[]').split(', ')]
    print(csv_files)

    dfs = [
        pd.read_csv(file, header=0, index_col=0)
        for file in csv_files
    ]
    df = pd.concat(dfs, axis=1)

    df.to_csv('merged.csv', index=True, header=True)
    """

}

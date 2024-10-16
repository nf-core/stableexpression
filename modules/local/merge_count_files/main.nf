process MERGE_COUNT_FILES {

    // label 'error_retry'
    debug true

    conda "${moduleDir}/environment.yml"

    input:
    val csv_files_str

    output:
    path 'all_counts.csv', emit: csv


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    from pathlib import Path

    # changing input value to list of paths
    csv_files = [Path(file) for file in "$csv_files_str".strip('[]').split(', ')]

    dfs = []
    for file in csv_files:

        # parsing dataframes
        df = pd.read_csv(file, header=0, index_col=0)

        # renaming columns, to avoid possible conflicts during concatenation (and especially subsequent modules)
        filename = file.stem.replace('_renamed', '').replace(',', '_').replace('.', '_')
        df.rename(columns={col: f'{filename}_{col}' for col in df.columns}, inplace=True)

        dfs.append(df)

    concat_df = pd.concat(dfs, axis=1)

    concat_df.to_csv('all_counts.csv', index=True, header=True)
    """

}

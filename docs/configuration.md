# Pipeline configuration

Although many parameters are directly exposed to the user via the CLI, setting cpu and memory requirements must be done via a configuration file.
By default, the max number of CPUs and memory the pipeline can use is defined in `conf/base.config`:

```
executor {
    cpus = 8
    memory = 24.GB
}
```

This was set quite low on purpose, in order to make it run easily on most data science laptops.

## Setting hard limits of CPU and memory

One can modify it by creating a custom config file:

```
executor {
    cpus = 32
    memory = 200.GB
}
```

then launch the pipeline using:

```bash
nextflow run nf-core/stableexpression \
   -profile <PROFILE> \
   -config <YOUR CONFIG FILE>
   ...
```

## Modifying resource allocation for processes

Let's say you have a laptop with only 4 CPUs and 12 GB of RAM. Running the pipeline may crash your computer because of a lack of memory.
To tell the pipeline to lower down its resource consumption, you can create a custom config file with:

```
executor {
    cpus = 4
    memory = 8.GB
}

withLabel:process_single {
    memory = { 2.GB * task.attempt }
}
withLabel:process_low {
    memory = { 2.GB + 1.GB * task.attempt }
}
withLabel:process_medium {
    memory = { 4.GB + 1.GB * task.attempt }
}
withLabel:process_high {
    memory = { 4.GB + 2.GB * task.attempt }
}
```

then launch the pipeline using:

```bash
nextflow run nf-core/stableexpression \
   -profile <PROFILE> \
   -config <YOUR CONFIG FILE>
   ...
```

> [!WARNING]
> Please keep in mind that if the total number of datasets (downloaded from public datasets or directly provided by the user) is too big for your computer, the pipeline will crash. Even if much effort was made to minimise the memory usage, some steps still require a certain amount of memory to run successfully.

## Running with Slurm

This is an example `launch_nf_core_stableexpression.sh` script to run the pipeline on an HPC cluster with Slurm:

```bash
#!/bin/bash

# set job name
#SBATCH --job-name=nf_run

# set output file for logs
#SBATCH --output=logs/nf_run_%j.log

# if your HPC cluster uses partitions, use a partition allowing for long runs
#SBATCH --partition=<LONG RUNTIME PARTITION>

#to get email notifications
#SBATCH --mail-user=<EMAIL>
#SBATCH --mail-type=END,FAIL

# set max memory available to the Nextflow main node
#SBATCH --mem 2GB

module load nextflow
# or load specific version: module load nextflow=25.10.04

# set location of apptainer cache directory
export NXF_APPTAINER_CACHEDIR=apptainer_cache

nextflow run nf-core/stableexpression \
    -latest \
    -c slurm.config \
    -profile apptainer \
    -resume \
    --params-file params.yaml
```

with `slurm.config`:

```
executor {
    name = 'slurm'
    queue = <FAST RUNTIME PARTITION> // if your HPC cluster uses partitions, use a partition including fast runs
    queueSize = 50 // see https://seqera.io/blog/5_tips_for_hpc_users/
    submitRateLimit = '10 sec' // see https://seqera.io/blog/5_tips_for_hpc_users/
    cpus = 64 // adjust to your needs
    memory = 400.GB // adjust to your needs
    time = 48.h // optional, only if you want to limit the runtime
}
```

and `params.yaml`:

```
species: <SPECIES>
outdir: <OUTDIR>
[+ OTHER PARAMETERS]
```

Run this script with `sbatch`:

```
sbatch launch_nf_core_stableexpression.sh
```

For checking the status of the run, we recommend tools like [slurmer](https://crates.io/crates/slurmer).

# Galaxy

## Setup build / testing environment

NB: You need conda installed (micromamba does not work, since the Galaxy installer looks for a venv / conda environment)

Create a new environment with python and planemo installed:

```
conda env create -f environment.yml -y
conda activate planemo
```

## Build tool XML file

The XML definition file is partially generated dynamically by:

- parsing nextflow_schema.json
- fetching latest version of Nextflow, Singularity and OpenJDK in Conda channels

However, you need to build a boilerplate file with things that cannot be directly interpreted from nextflow_schema.json, such as:

- path to selected output files
- tests
- specific conditions for the inputs

### Build boilerplate XML file (only once)

```
python build/build_boilerplate.py
```

The boilerplate XML file is generated at `galaxy/build/static/boilerplate.xml`.

### Customise boilerplate XML file

You must edit the boilerplate XML file to add your customisations:

- Mandatory (at least if your pipeline uses a samplesheet): modify file paths in the samplesheet
  Galaxy has its own path system, and you must retrieve dynamically the paths of the files provided, in order to modify them in the samplesheet
  "Running the pipeline"
  In this cas, add "&&" before "nextflow drop ..."

- modify outputs
- add tests

```
python build/build_custom.py
```

### Build XML file (at each release)

```
python build/build_tool.py
```

This script will fetch :

- all the parameters in your nextflow_schema.json
- the latest version of Nextflow, Singularity and OpenJDK in Conda channels.

Your tool is ready to be used!

## Test tool

### Launch local Galaxy server

You may want to have a first look at what your tool looks like in the Galaxy interface.
To launch a local instance of Galaxy with your tool already installed:

```
tool/serve.sh
```

You can test the behaviour of your tool by providing different inputs and check the corrsponding output.

### Linting and testing

To lint your tool:

```
test/lint.sh
```

To test your tool:

```
test/test.sh
```

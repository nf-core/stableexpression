# First time setup

>[!NOTE]
>The following instructions need only to be performed once when you want to initialise your Galaxy tool and your repository on the Galaxy Toolshed.

## Setup build / testing environment

Create a new environment with python and planemo installed:
```
micromamba env create -f environment.yml -y
micromamba activate planemo
```

## Initialise Galaxy tool boilerplate

The XML definition file is partially generated dynamically by:

- parsing nextflow_schema.json
- fetching latest version of Nextflow, Singularity and OpenJDK in Conda channels

However, you need to build a boilerplate file with things that cannot be directly interpreted from nextflow_schema.json, such as:

- path to selected output files
- tests
- specific conditions for the inputs

### Build template XML file

```
python build/build_boilerplate.py
```

The boilerplate XML file is generated at `galaxy/build/static/boilerplate.xml`.

### Customise template XML file

You must edit the boilerplate XML file to add your customisations:

- Mandatory (at least if your pipeline uses a samplesheet): modify file paths in the samplesheet
  Galaxy has its own path system, and you must retrieve dynamically the paths of the files provided, in order to modify them in the samplesheet
  "Running the pipeline"
  In this cas, add "&&" before "nextflow drop ..."

- modify outputs
- add tests

## Create repository on Toolshed

All necessary instructions are available in the [Galaxy Toolshed documentation](https://planemo.readthedocs.io/en/master/publishing.html).

For now, you just need to :
- [configure a shed account](https://planemo.readthedocs.io/en/master/publishing.html#configuring-a-shed-account)
- [create a new repository on the Toolshed](https://planemo.readthedocs.io/en/master/publishing.html#creating-a-repository)

Create a new folder for your tool and place the .shed.yml file in it:

```
mkdir -p tool_shed/tool
mv .shed.yml tool_shed
```

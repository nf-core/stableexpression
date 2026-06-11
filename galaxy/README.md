# Galaxy

The following instructions need to performed for each release

>[!TIP]
>For the first time setup of Galaxy for you Nextflow pipeline, see the [setup instructions](setup.md)

## Activate environment

>[!NOTE]
>If you're planemo environment is not set up, see the [setup instructions](setup.md)

Activate your planemo environment:
```
micromamba activate planemo
```

## At each release: build XML file

### Optional: modify static values in template file

If needed, you can:
- update the versions of core dependencies (Nextflow, Micromamba, OpenJDK)
- modify outputs
- modify tests

>[!NOTE]
>The versions of core dependencies (Nextflow, Micromamba, OpenJDK) are not updated automatically, although the code necessary for this is already implemented.
>For now, we want to keep control over the versions used, to avoid versions that may contain bugs.

### Update tool

Update the tool XML file:
```
python build/build_tool.py
```

This script will fetch :

- all the parameters in your nextflow_schema.json
- the latest version of Nextflow, Singularity and OpenJDK in Conda channels.

and modify the XML file located at `galaxy/tool_shed/tool/nf_core_stableexpression.xml`.

Your tool is ready to be used!

## Test tool

### Launch local Galaxy server

You may want to have a first look at what your tool looks like in the Galaxy interface.
To launch a local instance of Galaxy with your tool already installed:

```
./serve
```

You can test the behaviour of your tool by providing different inputs and check the corrsponding output.

### Linting and testing

To lint your tool:

```
./lint
```

>[!WARNING]
>The test script is not working for now... Planemo does not seem to find the input data for testing...
>For the moment, testing in a local webserver and linting using the provided script should be sufficient.

To test your tool:

```
./test
```

## Publishing to the Galaxy Toolshed


```
cd tool_shed
```

### Optional: test update on the test Toolshed

If you have already set up an account on the test Toolshed, you can test the update of your tool:
```
planemo shed_update --shed_target toolshed
```

### Official Galaxy Toolshed

```
planemo shed_update --shed_target toolshed
```

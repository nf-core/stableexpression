# Mimic runs of nf-test in Github runners using act

This folder contains all the necessary files to run `nf-test` tests using [act](https://nektosact.com/introduction.html).

## Install act

To install `act`, simply run:
```
curl --proto '=https' --tlsv1.2 -sSf https://raw.githubusercontent.com/nektos/act/master/install.sh | sudo bash
```

>[!NOTE]
>You might then have to place the act binary in a folder in your `$PATH`.

>[!IMPORTANT]
>`act` used `docker` under the hood. To install `docker`, see the [installation instructions](https://docs.docker.com/engine/install/).

## Setup tests to run

The `params.env` comprises all the necessary configuration to run the tests you need:
- profile(s)
- Nextflow version

## Run tests

You need to specify in `params.env` the profile(s) that will be used. All the other nf-test arguments must be provided as usual.

Example:
```.env
#params.env
NXF_VER=25.04.0
PROFILE=conda
```

```
# from the root folder of you repo
tests/act/run --tag <your tag> --debug --verbose
```

## Clean generated files
```
sudo rm -rf .nf-test 
```

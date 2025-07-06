# Galaxy

You need to create a virtual environment for running planemo. You cannot use conda for that (otherwise that would be too easy...)
See https://github.com/galaxyproject/galaxy/issues/20011

If venv and pip are not available on your system:

```bash
sudo apt install python3-venv python3-pip -y
```

NB: you may need to use pyenv to install a proper Python version

then create and setup your environment:

```bash
cd test
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

You need micromamba installed.

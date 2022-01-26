## Development

### Set up environment
```bash
python3.9 -m venv venv
source venv/bin/activate
pip install -U pip wheel setuptools
```

### Install dependencies

```bash
pip install pip-tools
pip-sync requirements-dev.txt
pre-commit install
```

### Run linters

```bash
# Pre-commit hooks
pre-commit run --all-files
```

### Update requirements

```bash
./dev-tools/rebuild-requirements.sh
```


# Batch processing


```bash
# Copy search results
python tools/copy-search-results.py


# Create Amber coordinate/topology files:
#   rst7, parm7 and corresponding pdb
python tools/prepare-structures.py


# Create MD protocols
python init.py


# Run MD protocols
python run_locally.py

```

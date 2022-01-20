## Development

### Set up environment
```bash
python3 -m venv venv
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

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


# Example of single structure `1wou` refinement using `LocalSlurmWorker` remote runner

1. One needs to download a PDB and structure factors CIF files from the RCSB database.

2. Next, convert CIF to MTZ using Phenix and expand it to P1 space group. Properly formatted results would be achieved with the following two commands:
```bash
phenix.cif_as_mtz 1WOU-sf.cif --ignore_bad_sigmas --merge --output_file_name=1wou-sf.mtz
phenix.reflection_file_converter --expand_to_p1 1wou-sf.mtz --write_mtz_amplitudes --mtz_root_label="FOBS" --label="FOBS" --generate_r_free_flags --non_anomalous --mtz 1wou.mtz
```

3. Create general use topology (and coordinates) file from PDB.
    * `source` AMBER, since we rely on its' Python libraries.
    * Put the PDB and P1-MTZ files into `data/input/1wou` directory of the cloned repo.
    * Setup Python virtual enviroment as discribed in this repo above. To avoid potential conflicts, use `amber.python`.
    * Add the repo's directory into `PYTHONPATH` variable (`export PYTHONPATH=$PYTHONPATH:$(pwd)`)
    * Run: `python tools/prepare-structures.py`. The results will be in `data/amber-topology/1wou`

4. Convert MTZ into a simple text format described in the AMBER manual and create x-ray specific topology file. The first part of this step is covered by the `write_sf_dat_file()` method of `init.py`. The second part of adjusting the general use topology file for x-ray purposes by adding x-ray and B-factors data is covered by the `prepare_xray_prmtop()` method. Thus, to prepare the simulation/refinement job run `python init.py` command. The results will be in `data/output/1wou`.

5. Finally, run the refinement job by executing `python run_remotely.py`! NB: don’t forget to adjust the paths that are `source`d and `export`ed. The results will be in `data/output/1wou`.

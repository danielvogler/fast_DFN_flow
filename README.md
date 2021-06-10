# fast_DFN_flow

### Needed
- overview test to ensure results do not change
- setup a conda environment
- translation to python3
- wrapper for graph-tool to invert dependency
- refactoring


### Conda environment:
```bash
conda env create -f environment.yml
```

This environment (Hobe2018etal) can then be activated using:

```bash
conda activate Hobe2018etal
```

To always have access to this kernel in jupyter notebooks you can add (not required, but useful):
```bash
python -m ipykernel install --user --name=indSeisHack
```

To test the addition of new packages to the environment.yml file :
```bash
conda env update --file environment.yml
```

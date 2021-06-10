# fast_DFN_flow

### Needed
- wrapper for graph-tool to invert dependency
- refactoring


### Dependencies
- graph-tool: [https://git.skewed.de/count0/graph-tool](https://git.skewed.de/count0/graph-tool)


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
python -m ipykernel install --user --name=Hobe2018etal
```

To test the addition of new packages to the environment.yml file :
```bash
conda env update --file environment.yml
```

### Unit testing:
The src/test folder includes tests to ensure the results do not change.
To check these tests, go into the src/test folder and use:
```bash
pytest -v
```
and
```bash
coverage report -m
```




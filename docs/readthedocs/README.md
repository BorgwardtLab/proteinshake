# Building docs

Install:

```
pip install sphinx sphinx-rtd-theme myst-nb torch torch_geometric
pip install git+https://github.com/BorgwardtLab/proteinshake.git
```
Then run to build:
```
make clean
make html
open build/html/index.html
```

Pushing to main branch automatically updates public docs.


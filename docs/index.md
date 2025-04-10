# Shiba

![SE](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/img/Shiba_icon.png?raw=true){ align=right width=40% }

A unified computational method for systematic identification of differential RNA splicing. **Shiba**/**scShiba** can quantify and identify differential splicing events (DSEs) from short-read bulk RNA-seq data and single-cell RNA-seq data. **Shiba** and **scShiba** are also implemented as [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows, **SnakeShiba** and **SnakeScShiba**, respectively.

See [CHANGELOG.md](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/CHANGELOG.md) for the latest updates.

!!! Important

    If you continue to encounter issues, please don't hesitate to [open an issue](https://github.com/Sika-Zheng-Lab/Shiba/issues) on GitHub. The community and developers are here to help!

## Contents

- [Installation](installation.md)
- Quick Start
    - [With bulk RNA-seq data](quickstart/diff_splicing_bulk.md)
    - [With single-cell RNA-seq data](quickstart/diff_splicing_sc.md)
- Output
    - [Shiba/SnakeShiba](output/shiba.md)
    - [scShiba/SnakeScShiba](output/scshiba.md)
- Usage
    - [shiba.py](usage/shiba.md)
    - [snakeshiba.smk](usage/snakeshiba.md)
    - [scshiba.py](usage/scshiba.md)
    - [snakescshiba.smk](usage/snakescshiba.md)

## Citation

Kubota N, Chen L, Zheng S. [Shiba: a versatile computational method for systematic identification of differential RNA splicing across platforms](https://academic.oup.com/nar/article/53/4/gkaf098/8042001). *Nucleic Acids Research*  53(4), 2025, gkaf098.

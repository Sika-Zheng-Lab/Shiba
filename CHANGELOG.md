# Change log

All notable changes to this Shiba project will be documented in this file.

## [v0.5.3] - 2025-03-05

### Added

- Print Shiba version in the log file.
- Generate a report for the version of pipeline, date of execution, and parameters used in the analysis.
- Deposit test datasets for Shiba and scShiba in Zenodo (https://zenodo.org/records/14976391) ([#40](https://github.com/Sika-Zheng-Lab/Shiba/issues/40#issuecomment-2660001877))

### Changed

- Update citation information in [README.md](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/README.md#citation).

### Fixed

- Fixed a bug in `bam2junc.py` due to wrong handling of duplicated juntions ([#40](https://github.com/Sika-Zheng-Lab/Shiba/issues/40)).
- Fixed a bug of misassignment of samples when there are >3 groups in an experiment table and calculating individual PSI values.

## [v0.5.2] - 2025-02-07

### Fixed

- Add pyYAML to the Dockerfile to fix the error in the execution of the main scripts and snakefiles.

## [v0.5.1] - 2025-01-21

### Added

- Added author information to README.md
  - Naoto Kubota ([0000-0003-0612-2300](https://orcid.org/0000-0003-0612-2300))
  - Liang Chen ([0000-0001-6164-4553](https://orcid.org/0000-0001-6164-4553))
  - Sika Zheng ([0000-0002-0573-4981](https://orcid.org/0000-0002-0573-4981))

## [v0.5.0] - 2025-01-15

> [!IMPORTANT]
> Main scripts were replaced with `shiba.py` and `scshiba.py` for Shiba and scShiba, respectively.
> The original scripts `Shiba` and `scShiba` were deprecated.

### Added

- Added `shiba.py` and `scshiba.py` as main scripts for Shiba and scShiba, respectively, and removed `Shiba` and `scShiba` to improve the usability of the software.
- Replaced bash scripts to Python scripts for better compatibility and maintainability.

### Changed

- Renamed snakefiles to `snakeshiba.smk` and `snakescshiba.smk` for SnakeShiba and SnakeScShiba, respectively.
- Improved logging messages for better readability.

### Fixed

- Optimized the code for better performance and readability.

## [v0.4.1] - 2024-10-16

### Added

### Changed

### Fixed

- Fix a bug in MSE detection when the number of MSEs is zero.

## [v0.4.0] - 2024-07-30

### Added

### Changed

- Add event type (e.g. `SE`, `FIVE`, `THREE`, `MXE`, `RI`, `MSE`, `AFE`, and `ALE`) to `pos_id` in events and results files.

### Fixed

- Fix a format of `pos_id` of `AFE` and `ALE` events to be consistent between different runs.

## [v0.3.1] - 2024-07-03

### Added

- Add GitHub actions for automatic tag creation and release, and Docker image build and push.

### Changed

### Fixed

## [v0.3] - 2024-07-02

### Added

- New features for identifying multiple skipped exons (MSE), alternative first exon (AFE), and alternative last exon (ALE) ([#4](https://github.com/Sika-Zheng-Lab/Shiba/issues/4))
  - Up to 5,000 tandemly skipped exons can be identified for MSE, which is reasonable for most datasets and organisms

### Changed

- Unannotated splicing events are now defined when at least one of the associated ***introns*** are not documented in the reference annotation file. Previously, unannotated splicing events were defined when at least one of the associated ***exons*** were not documented in the reference annotation file.
- Update [README.md](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/README.md) with new badges
- Update [MANUAL.md](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/doc/MANUAL.md) with new instructions for MSE, AFE, and ALE analyses

### Fixed

- Fixed minor bugs in `SnakeScShiba`

## [v0.2.2] - 2024-06-30

### Added

### Changed

### Fixed

- Fix `numpy` dependency in the Docker image to `1.26.4`

## [v0.2.1] - 2024-06-29

### Added

### Changed

### Fixed

- Fix the bug in `SnakeShiba` ([#7](https://github.com/Sika-Zheng-Lab/Shiba/issues/7))

## [v0.2] - 2024-06-17

### Added

- `pca.py`: Principal Component Analysis (PCA) module ([#1](https://github.com/Sika-Zheng-Lab/Shiba/issues/1))
  - Perform PCA for top 3,000 highly variable genes or splicing events
    - When PSI matrix has less than 6,000 rows without NaN values, drop rows with more than 50% NaN values followed by KNN imputation (n_neighbors=5), and then get top 3,000 highly variable splicing events
- This `CHANGELOG.md` file to keep track of changes

### Changed

- Generate PSI matrix (`PSI_matrix_sample.txt` and `PSI_matrix_group.txt`) when differential splicing analysis is performed
- Update the style of `summary.html`
- Update [MANUAL.md](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/doc/MANUAL.md) with new instructions
- Update [README.md](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/README.md) with citation information

### Fixed

## [v0.1] - 2024-05-17

### Added

- Initial release

### Changed

### Fixed

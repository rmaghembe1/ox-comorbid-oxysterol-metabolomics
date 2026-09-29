# OX-COMORBID: longitudinal oxysterol metabolomics

## Manuscript

**Calibration-constrained longitudinal oxysterol metabolomics reveals
coordinated remodeling of multivariate response architecture in
hypercholesterolaemia**

### Authors

- Reuben S. Maghembe
- Samweli Bahati
- Abdalah Makaranga
- Donath Damian

## Scope

This repository contains the computational reproducibility package for
a secondary analysis of previously generated and publicly available
digital metabolomics data.

The project contains no wet-laboratory experimentation, no collection or
handling of biological specimens, no organism or cell culture, no
genetic modification, and no experimental manipulation of biological
materials.

All analyses are restricted to existing digital analytical-response
data and published metadata.

## Analytical scope

The source panel contains five free oxysterol channels:

- 24S-hydroxycholesterol
- 25-hydroxycholesterol
- 27-hydroxycholesterol
- 7beta-hydroxycholesterol
- 7-ketocholesterol

The primary inferential framework uses calibration-robust quantities,
including paired directional topology, rank-space displacement,
rank-dependence rewiring and exact finite permutation spaces.

Response-scale analyses are retained as sensitivity analyses and are
not interpreted as participant-level concentration fold changes.

## Main result boundary

The analysis supports coordinated longitudinal remodeling of
multivariate oxysterol-response architecture.

The repository does **not** support claims that:

- the deposited corrected response values are absolute concentrations;
- the analysis establishes causal treatment effects;
- the statistical control reference is normative biological
  homeostasis;
- the five-analyte panel represents the entire oxysterolome.

## Repository structure

- `analysis/` - frozen computational analysis and result artifacts
- `figures/png/` - authoritative final PNG figures
- `figures/pdf/` - authoritative final PDF figures
- `publication_assets/tiff/` - authoritative high-resolution TIFF files
- `metadata/` - authorship, ORCID and contribution metadata
- `environment/` - runtime/version information
- `provenance/` - public package selection and integrity manifests

## Source data

The public source dataset is available at:

`10.17632/6cm8mxm5rm.1`

The raw source workbook is intentionally not redistributed in this
repository. See `SOURCE_DATA.md`.

## Release status

This is a local public-release candidate for version `v1.0.0`.

Public release remains **HOLD** until the repository licensing decision
is explicitly locked and the GitHub/Zenodo metadata are finalized.

## Licensing

This repository uses mixed licensing.

- Original executable analysis code: **MIT**
- Original project documentation and metadata: **CC BY 4.0**
- Source-derived analytical result surfaces and figures:
  **CC BY-NC 3.0**

See `LICENSE.md`, `LICENSE_MAP.tsv`, and
`SOURCE_LICENSE_AND_ATTRIBUTION.md`.

The original Mendeley workbook is not redistributed.

## Public release route

The planned public release route is:

1. GitHub repository;
2. immutable GitHub release `v1.0.0`;
3. Zenodo archival with explicit mixed-license metadata;
4. DOI verification before manuscript v3 is generated.

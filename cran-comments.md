# Resubmission

This is a resubmission. In this version I have:

* Reduced the size of the example datasets included in the package in order to reduce the size of the tarball.

* Added a reference to the publication presenting HIDECAN plots. It is currently in review,
so there is no DOI to cite yet. The package description will be updated upon publication.

# Test environments

- local Windows 10, R 4.2.2
- Winbuilder: x86_64-w64-mingw32 (64-bit) R version 4.2.2 (2022-10-31 ucrt)
- Winbuilder: x86_64-w64-mingw32 (64-bit) R Under development (unstable) (2023-01-25 r83685 ucrt)
- Rhub:
  - Windows Server 2022, R-devel, 64 bit
  - Fedora Linux, R-devel, clang gfortran
  - Ubuntu Linux 20.04 LTS, R-release, GCC
  - macOS 10.13.6 High Sierra, R-release, brew
  - debian-gcc-release: Debian Linux, R-release, GCC

# R CMD check results

There were no ERRORs, WARNINGs. 

There was 1 NOTE:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Olivia Angelin-Bonnet <olivia.angelin-bonnet@plantandfood.co.nz>'

New submission

Possibly misspelled words in DESCRIPTION:
  Angelin (11:19)
  GWAS (2:45, 9:18)
  HIDECAN (2:15, 8:24, 10:64)
  al (11:37)
  et (11:34)
  transcriptomics (9:27)
```

This is a new submission. Angelin is a family name; GWAS stands for Genome-Wide Association Study; HIDECAN is the name of the package; transcriptomics is correctly spelled; et al. was used for a reference.


# Downstream dependencies

There are currently no downstream dependencies for this package.

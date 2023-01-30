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
  GWAS (2:45, 9:18)
  HIDECAN (2:15, 8:24)
  transcriptomics (9:27)
```

This is a new submission. GWAS stands for Genome-Wide Association Study, HIDECAN is the name of the package, and transcriptomics is correctly spelled.

There was a NOTE only found on Windows Server 2022 R-devel, 64 bit:

```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```

As noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX and can likely be ignored.

There was a NOTE only found on macOS 10.13.6 High Sierra, R-release, brew:

```
* checking installed package size ... NOTE
  installed size is  5.2Mb
  sub-directories of 1Mb or more:
    extdata   1.0Mb
    R         3.1Mb
```

The extdata folder contains an R object that is made available to users as an example input data.

# Downstream dependencies

There are currently no downstream dependencies for this package.

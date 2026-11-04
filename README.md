# 🚨This repository is deprecated.

For the last QC implementations (2026), you can check the following studies:

- (Duran-Ferrer, Nat. Comms 2026, accepted)[https://github.com/Duran-FerrerM/CLL-ImmunoEpigenomics]. <br>
  QC using [SeSAMe](https://bioconductor.org/packages/3.23/bioc/vignettes/sesame/inst/doc/sesame.html) is avaiable [here](https://duran-ferrerm.github.io/CLL-ImmunoEpigenomics/revisions/revision-1/reviewer-1/Fig.1A.COO.submission.html)<br>

- (Gabbutt & Duran-Ferrer, Nature 2025)[https://github.com/Duran-FerrerM/evoflux]. <br>
  QC using using [minfi](https://bioconductor.org/packages/3.23/bioc/vignettes/minfi/inst/doc/minfi.html) is avaiable [here](https://duran-ferrerm.github.io/evoflux/QC_DNAme_arrays.v.4.1.html)
  
  **NOTE:** When using minfi, please mind the deafult SNP annotation db, as it may be outdated as previously reported:
  https://www.biostars.org/p/261603/ or https://www.biostars.org/p/9555147/.

  💡 In general, I recommend using SeSAME for its consistent update. However, note that the derivation of detection p-values is more stringent than the minfi procedure. You may consider this when merging several versions of Illumina DNA methylation arrays.

  ✉️ Do not hesitate to contact me for any doubt at maduran@recerca.clinic.cat.

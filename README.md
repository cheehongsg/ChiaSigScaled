# ChiaSigScaled-internal
ChIA-PET Interaction Caller (Internal)

## Overview
Jonas Paulsen published "[A statistical model of ChIA-PET data for accurate detection of chromatin 3D interactions](https://www.ncbi.nlm.nih.gov/pubmed/?term=25114054)" 
which provides a model parameterized by empirical data that improve the specificity of called significant ChIA-PET interactions.

However, the program could not scale to handle a single ChIA-PET library sequenced on a NextSeq. ChiaSigScaled is a heavily optimized implementation based on the original code by Paulsen.

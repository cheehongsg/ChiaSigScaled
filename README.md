# ChiaSigScaled
ChIA-PET Interaction Caller

## Overview
Jonas Paulsen published "[A statistical model of ChIA-PET data for accurate detection of chromatin 3D interactions](https://www.ncbi.nlm.nih.gov/pubmed/?term=25114054)" 
which provides a model parameterized by empirical data that improve the specificity of called significant ChIA-PET interactions.

The program, however, could not complete the processing of a ChIA-PET library output from a NextSeq. ChiaSigScaled is a heavily optimized implementation adapted from the [original code by Paulsen](http://folk.uio.no/jonaspau/chiasig/). A total of 44 optimization iterations have been completed to address various critical code hotspots. Optimization decisions have been made while maintaining the same results from Paulsen published paper dataset.


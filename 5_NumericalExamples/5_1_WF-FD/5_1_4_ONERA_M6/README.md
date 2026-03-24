# README

## Timings:

Runtime timings (comparison AVI vs. FD) are obtained on 32 threads on an AMD EPYC 9374F.

## Multirate Method:

The multirate method (`PairedExplicitRK4Multi`) and the pressure coefficient (`AnalysisSurfacePointwise`) are not implemented in a release of Trixi.jl. If you want to run them, reach out to me - I will point you to a branch where they are implemented.
On the same branch, the implementations are also much faster (almost factor of 3, also for the standalone schemes).

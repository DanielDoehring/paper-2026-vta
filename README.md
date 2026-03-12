# paper-2026-vta
Reproducibility Repository for the Paper "Volume Term Adaptivity for Discontinuous Galerkin Methods"


If you use the implementations provided here, please also cite this repository as
```bibtex
@misc{doehring2026VTA_ReproRepo,
  title={Reproducibility repository for "Volume Term Adaptivity for Discontinuous Galerkin Methods"},
  author={Doehring, Daniel and Jesse, Chan and Ranocha, Hendrik and Michael, Schlottke-Lakemper and Torrilhon, Manuel and Gregor, Gassner},
  year={2026},
  howpublished={\url{https://github.com/DanielDoehring/paper-2026-vta}},
  doi={https://doi.org/??}
}
```

## Abstract

We introduce the concept of volume term adaptivity for high–order discontinuous Galerkin (DG) schemes
solving time–dependent partial differential equations.
Termed *v*–adaptivity, we present a novel general approach that exchanges the discretization of the volume contribution
of the DG scheme at every Runge–Kutta stage based on suitable indicators.

Depending on whether robustness or efficiency is the main concern,
different adaptation strategies can be chosen.
Precisely, the weak form volume term discretization is used
instead of the entropy-conserving flux–differencing volume integral whenever the former produces more entropy than the latter, resulting in an entropy–stable scheme.
Conversely, if increasing the efficiency is the main objective,
the weak form volume integral may be employed as long as it does not increase
entropy beyond a certain threshold or cause instabilities.
Thus, depending on the choice of the indicator,
the *v*–adaptive DG scheme improves robustness, efficiency and approximation quality
compared to schemes with a uniform volume term discretization.

We thoroughly verify the accuracy, linear stability, and entropy–admissibility of the *v*-adaptive DG scheme to various compressible flow problems in two and three dimensions.

## Reproducing the results

### Installation

To download the code using `git`, use 

```bash
git clone git@github.com:DanielDoehring/paper-2026-vta.git
``` 

If you do not have git installed you can obtain a `.zip` and unpack it:
```bash
wget https://github.com/DanielDoehring/paper-2026-vta/archive/refs/heads/main.zip
unzip main.zip
mv paper-2026-vta-main/ paper-2026-vta
```

To instantiate the Julia environment, i.e., download all required packages, execute the following two commands:
```bash
cd paper-2026-vta/
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Note that the results are obtained using long-term-support (LTS) Julia 1.10.10, which is also set in the `Manifest.toml`.
Thus, you might need to install the [Julia 1.10.10 LTS release](https://julialang.org/downloads/manual-downloads/) first
and *replace* the `julia` calls from this README with
`/YOUR/PATH/TO/julia-1.10.10/bin/julia`

### Project initialization

If you installed Trixi.jl this way, you always have to start Julia with the `--project` flag set to your `paper-2026-vta` directory, e.g.,
```bash
julia --project=.
```
if already inside the `paper-2026-vta` directory.

If you do not execute from the `paper-2026-vta` directory, you have to call `julia` with
```bash
julia --project=/YOUR/PATH/TO/paper-2026-vta
```

### Running the code

The scripts for verification and applications are located in the `4_Verification`, and `5_NumericalExamples` directory, respectively.

To execute them provide the respective path:

```bash
julia --project=. ./4_Verification/4_1_Convergence/4_1_1_IsentropicVortex/elixir_euler_vortex_AVI.jl
```

For all cases in the `5_NumericalExamples` directory the solution has been computed using a specific number of threads.
To specify the number of threads the [`--threads` flag](https://docs.julialang.org/en/v1/manual/multi-threading/#Starting-Julia-with-multiple-threads) needs to be given, i.e., 
```bash
julia --project=. --threads 4 ./5_NumericalExamples/5_2_WF-FD-FV/5_2_1_RTI/elixir_euler_RTI_AMR_AVI.jl
```
The number of threads used for the examples are given in the `README.md` in `5_NumericalExamples`.

### Reproducing the figures

Plot scripts for selected plots (convergence, entropy, etc.) are provided in the respective directories.
The relevant datafiles are produced by the scripts.
For visualization of the simulation results, [`Trixi2Vtk](https://github.com/trixi-framework/Trixi2Vtk.jl) and [ParaView](https://www.paraview.org/) have been used.

## Authors

* [Daniel Doehring](https://www.acom.rwth-aachen.de/the-lab/team-people/name:daniel_doehring) (Corresponding Author)
* [Jesse Chan](https://sites.google.com/view/jessechan/home)
* [Hendrik Ranocha](https://ranocha.de/home#gsc.tab=0)
* [Michael Schlottke-Lakemper](https://lakemper.eu/)
* [Manuel Torrilhon](https://www.acom.rwth-aachen.de/the-lab/team-people/name:manuel_torrilhon)
* [Gregor Gassner](https://www.mi.uni-koeln.de/NumSim/gregor-gassner/)

## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
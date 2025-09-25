# Supernova Scalar Emission

This project computes the production and emission spectrum of CP-even scalars and the resulting luminosity of the scalar from a proto-neutron star (PNS) following a core-collapse supernova.  
The calculation uses precomputed PNS profiles (density, temperature, and chemical potentials) together with Monte Carlo phase-space integration.

---

## Requirements
- [Julia ≥ 1.9](https://julialang.org/)
- Julia packages:
  - `Cuba.jl`
  - `Interpolations.jl`
  - `DelimitedFiles.jl`
  - `Printf.jl`

Install dependencies from the Julia REPL:
```julia
using Pkg
Pkg.add(["Cuba", "Interpolations","Printf","DelimitedFiles"])

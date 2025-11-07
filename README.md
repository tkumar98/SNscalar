# Supernova Scalar Emission
This project computes the production and emission spectrum of **CP-even scalars** from a proto-neutron star (PNS) following a core-collapse supernova.  
Two production channels are included:

1. **Nucleon–nucleon bremsstrahlung**  
2. **Pion conversion processes**  

The file structure will be updated

---
## Requirements
- [Julia ≥ 1.9](https://julialang.org/)
- Julia packages:
  - `Cuba.jl`
  - `Interpolations.jl`
  - `DelimitedFiles.jl`
  - `Printf.jl`

Install dependencies:
```julia
using Pkg
Pkg.add(["Cuba", "Interpolations","Printf",""DelimitedFiles])
```
---

## PNS Profile Data
The code requires **time-dependent PNS profiles** (density, temperature, nucleon chemical potentials).  
These files are **not included in this repository**.

Publicly available profiles can be obtained from sources such as:
- [Garching Core-Collapse Supernova Archive]([https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/](https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/))
- Other published supernova simulation datasets.

---

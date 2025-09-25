# Supernova Scalar Emission
This project computes the production and emission spectrum of **CP-even scalars** from a proto-neutron star (PNS) following a core-collapse supernova.  
Two production channels are included:

1. **Nucleon–nucleon bremsstrahlung**  
2. **Pion conversion processes**  

Each channel has its own folder with code and outputs.

---

## Repository Structure
.
├── bremsstrahlung/ # Code + outputs for NN bremsstrahlung channel

├── pion_conversion/ # Code + outputs for pion conversion channel

├── README.md # Project documentation

---


Each folder contains:
- Julia code for scalar spectrum computation.
- Output spectra (`.dat` files) sorted by scalar energy.
- Example scripts for running with different parameters.

---

## Features
- Computes differential scalar emission spectrum  
  **d²N / (dt dE)** for each production channel.
- Flexible choice of scalar mass and mixing angle.
- Supports multiple time snapshots of PNS evolution.
- Outputs spectra in physical units (**s⁻¹ MeV⁻¹**).

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

Place the profile files inside a directory, for example:
SN_relevant/


---

## Usage

Run the main Julia script inside the desired channel folder:
```bash
cd bremsstrahlung/
julia spectrum.jl
```
---

Inside the script, specify desired times (in seconds after bounce):
```
user_times = [0.1, 1.0, 10.0]
```
---

Adjust Scalar Parameters
```
sin_theta = 1.0     # scalar mixing angle
ms = 1e-3           # scalar mass in GeV
```
---

Output

Each run produces spectra in text files:
```
scalar_spectrum_t{time}.dat
```

Columns:
```
Scalar energy (MeV)                Emission rate (s⁻¹ MeV⁻¹)
```

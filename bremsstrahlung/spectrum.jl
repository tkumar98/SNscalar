"""
SN_CPeven_Scalar_Spectrum.jl

This Julia script computes the spectrum of CP-even scalars emitted 
from a hot proto-neutron star (PNS) in a core-collapse supernova.

The code reads radial profiles of the PNS (temperature, density, 
chemical potentials, mean-field potentials, etc.) at different times 
from data files. It then computes the differential production rate 
of scalars via nucleon-nucleon and nucleon-pion interactions using 
Monte Carlo integration (Cuba VEGAS).

**Usage Notes for Public Release:**
- The PNS profile data files are proprietary and not included. 
  Users can generate similar data using publicly available PNS 
  simulation codes or get readymade ones from https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/.
- The user can now specify the simulation times of interest 
  instead of hardcoded indices.
- Units: GeV, km, cm, s, etc., are consistently used throughout.

Dependencies:
- Cuba.jl
- DelimitedFiles.jl
- Interpolations.jl
- Printf.jl
"""
using Cuba
using Printf
using DelimitedFiles
using Interpolations
using Base.Threads

# --------------------------- CONSTANTS ---------------------------

const R64 = Float64
const I64 = Int64
const Vec = Vector

# Unit conversions
const GeV_cm_inv = 5.06e13    # GeV -> cm^-1
const cm_inv_GeV = 1.0 / GeV_cm_inv
const cm_GeV_inv = cm_inv_GeV
const muG_GeV = 1.95e-26     # microGauss -> GeV^2
const km_Mpc = 1.0 / 3.09e19 # km -> Mpc
const kpc_cm = 3.086e21      # kpc -> cm
const cm_TeV_inv = 1.e3 / cm_inv_GeV
const gm_GeV = 5.62e23
const G_GeV = 1.95e-20
const keV_GeV = 1.e-6
const GeV_keV = 1.e6
const eV_GeV = 1.e-9
const GeV_erg = 1.6e-3
const GeV_s_inv = 1.52e24
const K_GeV = 8.62e-14
const erg_GeV = 6.24e2
const s_inv_GeV = 6.58e-25
const s_GeV_inv = 1.52e24
const yr_s = 3.16e7
const km_cm = 1.e5
const fm_cm = 1.e-13

# Nuclear and particle physics constants
const mnucl = 0.938919      # Nucleon mass (GeV)
const mn = 0.939
const mp = mn - 1.293e-3
const yhNN = 1.e-3
const vEW = 246.0           # Electroweak vev (GeV)
const mpi = 134.9768 / 1.e3
const mpiplus = 139.57039 / 1.e3
const gA = 1.28
const fpi = 92.4e-3
const nB = 1.2e38            # Baryon number density (cm^-3)
const TSN0 = 30.e-3          # Typical SN temperature (GeV)

# --------------------------- HELPER FUNCTIONS ---------------------------

"""Pion propagator function f(psq) = 1/(psq + m_pi^2)"""
f(psq) = 1.0 / (psq + mpi^2)

"""Squared pion propagator"""
g(psq) = 1.0 / (psq + mpi^2)^2

"""Charged pion propagator"""
fprime(psq) = 1.0 / (psq + mpiplus^2)

"""Squared charged pion propagator"""
gprime(psq) = 1.0 / (psq + mpiplus^2)^2

"""Fermi-Dirac distribution"""
f_fermi(x, T) = 1.0 / (exp(x / T) + 1)

# --------------------------- LOAD PNS PROFILES ---------------------------

# Path to PNS profile data files (user should replace with actual path)
path_to_profiles = "SN_relevant"

# Collect all .dat files
files = filter(f -> endswith(f, ".dat"), readdir(path_to_profiles))

# Extract simulation times from filenames (pattern: 0p12 -> 0.12)
sim_time_points = Float64[]
for file in files
    m = match(r"\d+p\d+", file)
    if m !== nothing
        value = parse(Float64, replace(m.match, "p" => "."))
        push!(sim_time_points, value)
    end
end

# --------------------------- USER INPUT ---------------------------

"""
Specify the times of interest (in seconds) for which to compute the 
scalar spectra. The code will automatically match the nearest available 
simulation time point from the PNS data files.
"""
println("Available simulation times (s): ", sim_time_points)
println("Enter desired times as an array of floats, e.g., [0.1, 1.2, 10.0]:")
user_times = parse.(Float64, split(readline(), ","))

# Find indices of closest available time points
using Statistics
t_index_arr = [findmin(abs.(sim_time_points .- t))[2] for t in user_times]

# --------------------------- MAIN LOOP OVER TIMES ---------------------------

@threads for index in t_index_arr
    t = sim_time_points[index]
    time_str = replace(@sprintf("%.8f", t), "." => "p")
    input_filename = time_str * ".dat"
    input_filepath = joinpath(path_to_profiles, input_filename)

    # Read profile data (columns: radius, T, rho, alpha, mu_e, mu_mu, mu_n, mu_p, mu_pi, U_n, U_p, mn*, mp*)
    data = readdlm(input_filepath, skipstart=1)
    radii = data[:, 1]

    # Interpolate PNS profiles for continuous radius values
    temp_R     = LinearInterpolation(radii, data[:, 2], extrapolation_bc=Flat())      # GeV
    rho_R      = LinearInterpolation(radii, data[:, 3] / 1.e14, extrapolation_bc=Flat()) # 10^14 g/cm^3
    alpha_R    = LinearInterpolation(radii, data[:, 4], extrapolation_bc=Flat())      # dimensionless
    mu_e_R     = LinearInterpolation(radii, data[:, 5], extrapolation_bc=Flat())      # GeV
    mu_mu_R    = LinearInterpolation(radii, data[:, 6], extrapolation_bc=Flat())      # GeV
    mu_n_R     = LinearInterpolation(radii, data[:, 7], extrapolation_bc=Flat())      # GeV
    mu_p_R     = LinearInterpolation(radii, data[:, 8], extrapolation_bc=Flat())      # GeV
    mu_pi_R    = LinearInterpolation(radii, data[:, 9], extrapolation_bc=Flat())      # GeV
    U_n_R      = LinearInterpolation(radii, data[:, 10], extrapolation_bc=Flat())     # GeV
    U_p_R      = LinearInterpolation(radii, data[:, 11], extrapolation_bc=Flat())     # GeV
    mn_star_R  = LinearInterpolation(radii, data[:, 12], extrapolation_bc=Flat())     # GeV
    mp_star_R  = LinearInterpolation(radii, data[:, 13], extrapolation_bc=Flat())     # GeV

    output_filename = "spectrum_" * time_str * ".dat"
    output_file = open(output_filename, "w")

    # --------------------------- SCALAR PRODUCTION INTEGRAND ---------------------------

    """
    S_integrand(x, func, params)

    Computes the integrand for the differential scalar production rate
    using nucleon-nucleon and nucleon-pion interactions.

    Inputs:
    - x: 6-dimensional integration variables (radius, momenta, angles)
    - func: array of length 1 to store output integrand
    - params: [sin_theta, ms, Es_inf] with:
        sin_theta: scalar mixing angle
        ms: scalar mass (GeV)
        Es_inf: scalar energy at infinity (GeV)
    """
    function S_integrand(x, func, params)
        sin_theta, ms, Es_inf = params[1], params[2], params[3]

        # Energetically forbidden
        if Es_inf < ms
            func[1] = 0.0
            return
        end

        # Map integration variable to physical radius
        xrmin, xrmax = 0.11666667, 20.0   # km
        dxr = xrmax - xrmin
        xr = xrmin + x[1] * dxr
        R = xr

        # Interpolate PNS profiles at radius R
        TSN = temp_R(R)
        alpha = alpha_R(R)
        mu_n = mu_n_R(R)
        mu_p = mu_p_R(R)
        U_n = U_n_R(R)
        U_p = U_p_R(R)
        mn_star = mn_star_R(R)
        mp_star = mp_star_R(R)

        # Scale integration variables
        dx1, dx2, dx3 = 10.0, 10.0, 10.0
        dcos_theta13, dcos_theta12 = 2.0, 2.0

        x1, x2, x3 = dx1 * x[2], dx2 * x[3], dx3 * x[4]
        cos_theta13 = -1.0 + dcos_theta13 * x[5]
        cos_theta12 = -1.0 + dcos_theta12 * x[6]

        mnucl_star = 0.5 * (mn_star + mp_star)
        E5 = Es_inf / alpha

        r, q, y = mnucl / TSN, ms / TSN, E5 / TSN
        sin_theta13, sin_theta12 = sqrt(1.0 - cos_theta13^2), sqrt(1.0 - cos_theta12^2)
        p1, p2, p3 = x1 * sqrt(mnucl * TSN), x2 * sqrt(mnucl * TSN), x3 * sqrt(mnucl * TSN)
        E1, E2, E3 = mnucl + p1^2 / 2.0 / mnucl, mnucl + p2^2 / 2.0 / mnucl, mnucl + p3^2 / 2.0 / mnucl

        # Compute angle beta
        cos_beta = -(
            (2.0 * mnucl^2 + E5^2 - 2.0 * mnucl * E5 +
             2.0 * (E1*E2 - p1*p2*cos_theta12) -
             2.0 * (E1*E3 - p1*p3*cos_theta13) -
             2.0 * (E2*E3 - p2*p3*cos_theta12*cos_theta13)
            ) / (2.0 * p2 * p3 * sin_theta12 * sin_theta13)
        )

        if cos_beta^2 > 1.0
            func[1] = 0.0
            return
        end

        E4 = E1 + E2 - E3 - E5

        # Kinematic factors for integrand
        abar = -1.0 * x2^2 * (x1^2 + x3^2 - 2.0 * x1 * x3 * cos_theta13)
        bbar = x2 * (x1 - x3 * cos_theta13) * (y^2 / r - 2.0*y - 2.0*x3^2 + 2.0*x1*x3*cos_theta13)
        cbar = x2^2 * x3^2 * sin_theta13^2 - (y^2/2.0/r - y - x3^2 + x1*x3*cos_theta13)^2
        Abar = abar * cos_theta12^2 + bbar * cos_theta12 + cbar

        if Abar <= 0.0 || E4 < mnucl
            func[1] = 0.0
            return
        end

        # ------------------ Matrix elements ------------------
        """
        Modify it according to the specific CP-even scalar model.
        Here we have adopted the Lagrangian provided in https://arxiv.org/abs/2005.00490
        """
        ksq = p1^2 + p3^2 - 2.0*p1*p3*cos_theta13
        lsq = p2^2 + p3^2 - 2.0*p2*p3*(cos_theta12*cos_theta13 + sin_theta12*sin_theta13*cos_beta)
        kdotl = p2*p3*(cos_theta12*cos_theta13 + sin_theta12*sin_theta13*cos_beta) - p3^2 - p1*p2*cos_theta12 + p1*p3*cos_theta13

        g1 = yhNN*sin_theta
        g2 = (2.0/9.0/vEW) * (ms^2 + 5.5*mpi^2) * sin_theta
        g4 = gA / (2.0*fpi)

        delta = ms^2 / (2.0*mnucl*E5)
        beta1 = 16.0*mnucl^2*g1*g4^2*delta / E5
        beta2 = 8.0*mnucl^2*g2*g4^2

        MNNsq = (1/4) * (
            (beta1*f(ksq) + beta2*g(ksq))^2 * 4.0 * ksq^2 +
            (beta1*f(lsq) + beta2*g(lsq))^2 * 4.0 * lsq^2 -
            (beta1*f(ksq) + beta2*g(ksq))*(beta1*f(lsq) + beta2*g(lsq))*(-4.0)*(ksq*lsq - 2.0*kdotl^2)
        )

        Mnpsq = (
            (beta1*f(ksq) + beta2*g(ksq))^2 * 4.0 * ksq^2 +
            (2.0*beta1*fprime(lsq) + beta2*gprime(lsq))^2 * 4.0 * lsq^2 +
            (2.0*beta1*fprime(lsq) + beta2*gprime(lsq))*(beta1*f(ksq) + beta2*g(ksq))*(-4.0)*(ksq*lsq - 2.0*kdotl^2)
        )

        ratio_mnucl = mnucl / mnucl_star

        # ------------------ Distribution functions ------------------
        eta_n = mu_n - mnucl - U_n
        eta_p = mu_p - mnucl - U_p

        Fnn = f_fermi(0.5*TSN*ratio_mnucl*x1^2 - eta_n, TSN) *
              f_fermi(0.5*TSN*ratio_mnucl*x2^2 - eta_n, TSN) *
              (1.0 - f_fermi(0.5*TSN*ratio_mnucl*x3^2 - eta_n, TSN)) *
              (1.0 - f_fermi(0.5*TSN*ratio_mnucl*(x1^2+x2^2-x3^2)-E5-eta_n, TSN))

        Fpp = f_fermi(0.5*TSN*ratio_mnucl*x1^2 - eta_p, TSN) *
              f_fermi(0.5*TSN*ratio_mnucl*x2^2 - eta_p, TSN) *
              (1.0 - f_fermi(0.5*TSN*ratio_mnucl*x3^2 - eta_p, TSN)) *
              (1.0 - f_fermi(0.5*TSN*ratio_mnucl*(x1^2+x2^2-x3^2)-E5-eta_p, TSN))

        Fnp = f_fermi(0.5*TSN*ratio_mnucl*x1^2 - eta_n, TSN) *
              f_fermi(0.5*TSN*ratio_mnucl*x2^2 - eta_p, TSN) *
              (1.0 - f_fermi(0.5*TSN*ratio_mnucl*x3^2 - eta_n, TSN)) *
              (1.0 - f_fermi(0.5*TSN*ratio_mnucl*(x1^2+x2^2-x3^2)-E5-eta_p, TSN))

        # ------------------ Integrand ------------------
        integ = (1.0 / alpha) * (TSN / TSN0)^4.5 * sqrt(y^2 - q^2) *
                xr^2 * x1^2 * x2^2 * x3^2 *
                (MNNsq*Fnn + MNNsq*Fpp + Mnpsq*Fnp) *
                (1.0 / sqrt(Abar)) * dx1 * dx2 * dx3 * dcos_theta12 * dcos_theta13 * dxr

        func[1] = integ  # GeV^-2
        return
    end

    # --------------------------- DIFFERENTIAL SPECTRUM FUNCTION ---------------------------
    function d2NdtdE(sin_theta, ms, Es_inf)
        """
        Computes the differential production rate of scalars (d^2N/dt dE)
        using VEGAS Monte Carlo integration.

        Inputs:
        - sin_theta: scalar mixing angle
        - ms: scalar mass (GeV)
        - Es_inf: scalar energy at infinity (GeV)

        Returns:
        - S: differential rate (dimensionless, can be converted to s^-1 MeV^-1)
        """
        params = [sin_theta, ms, Es_inf]

        # Prefactors: temperature, nucleon mass, volume
        prefact = 0.5 * TSN0^4.5 * mnucl^0.5 / (2*pi)^7
        volume_prefact = (1.0 * km_cm * cm_GeV_inv)^3

        # VEGAS integration parameters
        CUBA_minevals, CUBA_maxevals = 1_000_000, 100_000_000
        CUBA_nstart, CUBA_nincrease = 1_000_000, 1_000_000
        CUBA_rtol = 1.0e-3

        res = Cuba.vegas(
            S_integrand, 6, 1, userdata=params,
            rtol=CUBA_rtol,
            minevals=CUBA_minevals, maxevals=CUBA_maxevals,
            nstart=CUBA_nstart, nincrease=CUBA_nincrease
        )

        S = prefact * volume_prefact * res.integral[1]
        return S
    end

    # --------------------------- ENERGY LOOP ---------------------------
    # Compute the scalar emission spectrum (d^2N/dt dE) over a range of energies.

    sin_theta = 1.0      # Scalar mixing angle
    ms = 1.e-3           # Scalar mass in GeV

    # Log-spaced energies from 1e-3 to 1 GeV (50 points)
    Es_arr = 10.0 .^ LinRange(-3.0, 0.0, 50)

    for Es_inf in Es_arr
        # Compute the differential emission rate for this energy
        spec_val = d2NdtdE(sin_theta, ms, Es_inf)

        # Convert to physical units: s^-1 MeV^-1
        spec_dim = spec_val * GeV_s_inv * 1.e-3

        # Write energy (MeV) and spectrum to output file
        @printf(output_file, "%.6e  %.6e\n", Es_inf*1.e3, spec_dim)
        flush(output_file)  # ensure data is written immediately
    end

    close(output_file)

    # Sort output by energy
    run(`sort -g -k1,1 $output_filename -o $output_filename`)
end


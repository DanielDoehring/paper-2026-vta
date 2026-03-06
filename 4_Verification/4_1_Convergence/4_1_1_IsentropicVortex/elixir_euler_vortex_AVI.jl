using OrdinaryDiffEqLowStorageRK
using Trixi

# Ratio of specific heats
gamma = 1.4
equations = CompressibleEulerEquations2D(gamma)

surface_flux = flux_hllc

volume_flux = flux_ranocha
polydeg = 3
basis = LobattoLegendreBasis(polydeg)

volume_integral_weakform = VolumeIntegralWeakForm()
volume_integral_fluxdiff = VolumeIntegralFluxDifferencing(volume_flux)

# This indicator compares the entropy production of the weak form to the 
# entropy-conserving flux-differencing volume integral.
# If the entropy production of the weak form is lower than that of the
# flux-differencing form, we use the flux-differencing form to stabilize the solution.
indicator = IndicatorEntropyChange(maximum_entropy_increase = 0.0)

# Adaptive volume integral using the entropy production comparison indicator to perform the 
# stabilized/EC volume integral when needed and keeping the weak form if it is more diffusive.
volume_integral = VolumeIntegralAdaptive(volume_integral_default = volume_integral_weakform,
                                         volume_integral_stabilized = volume_integral_fluxdiff,
                                         indicator = indicator)

solver = DGSEM(basis, surface_flux, volume_integral)

"""
    initial_condition_isentropic_vortex(x, t, equations::CompressibleEulerEquations2D)

The classical isentropic vortex test case as presented in Section 5.1 of

- Brian Vermeire (2019).
  Paired Explicit Runge-Kutta Schemes for Stiff Systems of Equations
  [DOI:10.1016/j.jcp.2019.05.014](https://doi.org/10.1016/j.jcp.2019.05.014)
  https://spectrum.library.concordia.ca/id/eprint/985444/1/Paired-explicit-Runge-Kutta-schemes-for-stiff-sy_2019_Journal-of-Computation.pdf
"""
function initial_condition_isentropic_vortex(x, t, equations::CompressibleEulerEquations2D)
    # Evaluate error after full domain traversion
    if t == t_end()
        t = zero(t)
    end

    RealT = eltype(x)
    # Initial center of the vortex
    inicenter = SVector(0, 0)
    # Strength of the vortex
    S = convert(RealT, 13.5)
    # Radius of vortex
    R = convert(RealT, 1.5)
    # Free-stream Mach 
    M = convert(RealT, 0.4)
    # Base flow
    v1 = 1
    v2 = 1
    vel = SVector(v1, v2)

    center = inicenter + vel * t # Advection of center
    center = x - center          # Distance to centerpoint
    center = SVector(center[2], -center[1])
    r2 = center[1]^2 + center[2]^2

    f = (1 - r2) / (2 * R^2)

    rho = (1 - (S * M / convert(RealT, pi))^2 * (gamma - 1) * exp(2 * f) / 8)^(1 /
                                                                               (gamma - 1))

    du = S / (2 * convert(RealT, pi) * R) * exp(f) # Vel. perturbation
    vel = vel + du * center
    v1, v2 = vel

    p = rho^gamma / (gamma * M^2)

    prim = SVector(rho, v1, v2, p)
    return prim2cons(prim, equations)
end
initial_condition = initial_condition_isentropic_vortex

edge_length() = 20.0

N_passes() = 1
t_end() = edge_length() * N_passes()
tspan = (0.0, t_end())

coordinates_min = (-edge_length() / 2, -edge_length() / 2)
coordinates_max = (edge_length() / 2, edge_length() / 2)

mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 4,
                n_cells_max = 100_000,
                periodicity = true)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_condition_periodic)

ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (entropy,))

stepsize_callback = StepsizeCallback(cfl = 2.0)

callbacks = CallbackSet(summary_callback,
                        stepsize_callback,
                        analysis_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false, thread = Trixi.True());
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks);

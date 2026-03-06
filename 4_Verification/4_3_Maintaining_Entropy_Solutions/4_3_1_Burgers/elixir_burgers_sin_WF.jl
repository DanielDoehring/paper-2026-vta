using OrdinaryDiffEqSSPRK
using Trixi

###############################################################################
# semidiscretization of the (inviscid) Burgers' equation

equations = InviscidBurgersEquation1D()

basis = LobattoLegendreBasis(3)

surface_flux = flux_godunov
solver = DGSEM(basis, surface_flux, VolumeIntegralWeakForm())

coordinate_min = 0.0
coordinate_max = 1.0

# Make sure to turn periodicity explicitly off as special boundary conditions are specified
mesh = TreeMesh(coordinate_min, coordinate_max,
                initial_refinement_level = 6,
                n_cells_max = 10_000,
                periodicity = true)

function initial_condition_shock(x, t, equation::InviscidBurgersEquation1D)
    scalar = sinpi(2 * x[1]) + 0.5

    return SVector(scalar)
end


initial_condition = initial_condition_shock

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_condition_periodic)

###############################################################################
# ODE solvers, callbacks etc.

t_discontinuous = 1/(2 * pi) # time when shock forms

t_end = 0.25
#t_end = 10.0 # test long time stability
tspan = (0.0, t_end)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

stepsize_callback = StepsizeCallback(cfl = 1.2)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, SSPRK54(); dt = 1.0,
            ode_default_options()..., callback = callbacks);

using Plots

# Plot for weak form
plot(sol,
     guidefont = font("Computer Modern", 16), tickfont = font("Computer Modern", 14),
     titlefont = font("Computer Modern", 18), legendfont = font("Computer Modern", 16),
     labelfont = font("Computer Modern", 14),
     linewidth = 2, color = RGB(246/256, 169/256, 0),
     label = "Weak Form",
     #title = "Shock Formation Burgers'",
     title = "",
     legend = :bottomleft,
     yticks = [-20, -10, 0, 10, 20], ylim = (-21, 21),
     xlim = (0, 1.02),
     ylabel = "\$u\$",
     dpi = 600)


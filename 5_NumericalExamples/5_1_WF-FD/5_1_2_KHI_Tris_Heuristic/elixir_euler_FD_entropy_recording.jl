using OrdinaryDiffEqLowStorageRK
using Trixi

volume_integral_fluxdiff = VolumeIntegralFluxDifferencing(flux_ranocha)

dg = DGMulti(polydeg = 3,
             element_type = Tri(), # `Tri()` makes flux differencing really(!) expensive
             approximation_type = Polynomial(),
             surface_integral = SurfaceIntegralWeakForm(flux_hllc),
             volume_integral = volume_integral_fluxdiff)

equations = CompressibleEulerEquations2D(1.4)

"""
    initial_condition_kelvin_helmholtz_instability(x, t, equations::CompressibleEulerEquations2D)

A version of the classical Kelvin-Helmholtz instability based on
- Andrés M. Rueda-Ramírez, Gregor J. Gassner (2021)
  A Subcell Finite Volume Positivity-Preserving Limiter for DGSEM Discretizations
  of the Euler Equations
  [arXiv: 2102.06017](https://arxiv.org/abs/2102.06017)
"""
function initial_condition_kelvin_helmholtz_instability(x, t,
                                                        equations::CompressibleEulerEquations2D)
    # change discontinuity to tanh
    # typical resolution 128^2, 256^2
    # domain size is [-1,+1]^2
    slope = 15
    amplitude = 0.02
    B = tanh(slope * x[2] + 7.5) - tanh(slope * x[2] - 7.5)
    rho = 0.5 + 0.75 * B
    v1 = 0.5 * (B - 1)
    v2 = 0.1 * sin(2 * pi * x[1])
    p = 1.0
    return prim2cons(SVector(rho, v1, v2, p), equations)
end
initial_condition = initial_condition_kelvin_helmholtz_instability

cells_per_dimension = (32, 32)
mesh = DGMultiMesh(dg, cells_per_dimension; periodicity = true)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, dg;
                                    boundary_conditions = boundary_condition_periodic)

tspan = (0.0, 3.4) # FD still stable until this time

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

stepsize_callback = StepsizeCallback(cfl = 1.0)

analysis_interval = 10
analysis_callback = AnalysisCallback(semi, interval = analysis_interval, uEltype = real(dg),
                                     analysis_errors = Symbol[],
                                     extra_analysis_integrals = (entropy,),
                                     output_directory = "5_NumericalExamples/5_1_WF-FD/5_1_2_KHI_Tris_Heuristic/",
                                     analysis_filename = "analysis_FD.dat",
                                     save_analysis = true)

callbacks = CallbackSet(summary_callback,
                        stepsize_callback,
                        analysis_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0, ode_default_options()..., callback = callbacks);

using Trixi
using LinearAlgebra: eigvals
using OrdinaryDiffEqLowStorageRK

###############################################################################

equations = CompressibleEulerEquations2D(1.4)

# Setup from:
#=
- Gregor J. Gassner, Magnus Svärd, Florian J. Hindenlang (2020)
  Stability issues of entropy-stable and/or split-form high-order schemes
  [arXiv: 2007.09026](https://arxiv.org/abs/2007.09026)
=#

surface_flux = flux_lax_friedrichs

polydeg = 5
basis = LobattoLegendreBasis(polydeg)

volume_integral_weakform = VolumeIntegralWeakForm()

volume_integral = volume_integral_weakform
solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (-1.0, -1.0)
coordinates_max = (1.0, 1.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 2, # 4 x 4
                n_cells_max = 100_000,
                periodicity = true)

initial_condition = initial_condition_density_wave
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_condition_periodic)

J = jacobian_ad_forward(semi)
Eigenvalues = eigvals(J)
println("Max real value: ", maximum(real.(Eigenvalues)))


EigValFile = "./4_Verification/4_2_LinearStability/Eigenvalues_t0_WF.txt"
ofstream = open(EigValFile, "w")
for i in eachindex(Eigenvalues)
  realstring = string(real(Eigenvalues[i]))
  write(ofstream, realstring)

  write(ofstream, "+")

  imstring = string(imag(Eigenvalues[i]))
  write(ofstream, imstring)
  write(ofstream, "i") # Cpp uses "I" for the imaginary unit
  if i != length(Eigenvalues)
    write(ofstream, "\n")
  end
end
close(ofstream)


tspan = (0.0, 5.0) # 5 or 5000
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 5
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                      save_analysis = true,
                                      extra_analysis_integrals = (entropy,),
                                      analysis_errors = Symbol[],
                                      output_directory="./4_Verification/4_2_LinearStability/",
                                      analysis_filename = "analysis_WF.dat")

alive_callback = AliveCallback(alive_interval = 100)

stepsize_callback = StepsizeCallback(cfl = 0.9) # 0.05 in Paper mentioned above used

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 5e-5,
            ode_default_options()..., callback = callbacks);


using Trixi
using ForwardDiff
using LinearAlgebra: eigvals
using OrdinaryDiffEqLowStorageRK

# Required for AD jacobian

@inline function volume_integral_kernel!(du::AbstractArray{<:ForwardDiff.Dual}, u, 
                                         element, mesh,
                                         have_nonconservative_terms, equations,
                                         volume_integral::VolumeIntegralAdaptive{<:IndicatorEntropyChange},
                                         dg::DGSEM, cache)
    @unpack volume_integral_default, volume_integral_stabilized, indicator = volume_integral
    @unpack maximum_entropy_increase = indicator

    volume_integral_kernel!(du, u, element, mesh,
                            have_nonconservative_terms, equations,
                            volume_integral_default, dg, cache)

    # Compute entropy production of the default volume integral.
    # Minus sign because of the flipped sign of the volume term in the DG RHS.
    # No scaling by inverse Jacobian here, as there is no Jacobian multiplication
    # in `integrate_reference_element`.
    dS_default = -entropy_change_reference_element(du, u, element,
                                                   mesh, equations, dg, cache)

    # Compute true entropy change given by surface integral of the entropy potential
    dS_true = surface_integral_reference_element(entropy_potential, u, element,
                                                 mesh, equations, dg, cache)

    entropy_change = dS_default - dS_true

    # NOTE: For AD: Need `value` here
    if entropy_change.value > maximum_entropy_increase.value # Recompute using EC FD volume integral
        # Reset default volume integral contribution.
        # Note that this assumes that the volume terms are computed first,
        # before any surface terms are added.
        du[.., element] .= zero(eltype(du))

        volume_integral_kernel!(du, u, element, mesh,
                                have_nonconservative_terms, equations,
                                volume_integral_stabilized, dg, cache)
    end

    return nothing
end

###############################################################################

equations = CompressibleEulerEquations2D(1.4)

# Setup from:
#=
- Gregor J. Gassner, Magnus Svärd, Florian J. Hindenlang (2020)
  Stability issues of entropy-stable and/or split-form high-order schemes
  [arXiv: 2007.09026](https://arxiv.org/abs/2007.09026)
=#

surface_flux = flux_lax_friedrichs

volume_flux = flux_chandrashekar

polydeg = 5
basis = LobattoLegendreBasis(polydeg)

volume_integral_weakform = VolumeIntegralWeakForm()
volume_integral_fluxdiff = VolumeIntegralFluxDifferencing(volume_flux)

# This indicator compares the entropy production of the weak form to the
# true entropy evolution in that cell.
# If the weak form dissipates more entropy than the true evolution
# the indicator renders this admissible. Otherwise, the more stable
# volume integral is to be used.
indicator = IndicatorEntropyChange(maximum_entropy_increase = 0.0)

# Adaptive volume integral using the entropy production comparison indicator to perform the
# stabilized/EC volume integral when needed and keeping the weak form if it is more diffusive.
volume_integral_adaptive = VolumeIntegralAdaptive(volume_integral_default = volume_integral_weakform,
                                                  volume_integral_stabilized = volume_integral_fluxdiff,
                                                  indicator = indicator)

volume_integral = volume_integral_adaptive

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


EigValFile = "./4_Verification/4_2_LinearStability/Eigenvalues_t0_AVI.txt"
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
                                      analysis_filename = "analysis_AVI.dat")

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


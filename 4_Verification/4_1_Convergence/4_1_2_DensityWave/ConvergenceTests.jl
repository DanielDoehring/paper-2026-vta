using Trixi

cd(dirname(@__FILE__))

convergence_test("elixir_euler_density_wave_AVI.jl", 7)
convergence_test("elixir_euler_density_wave_FD.jl", 7)
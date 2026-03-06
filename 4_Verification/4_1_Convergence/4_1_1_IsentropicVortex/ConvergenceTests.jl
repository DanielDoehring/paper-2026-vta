using Trixi

cd(dirname(@__FILE__))

convergence_test("elixir_euler_vortex_AVI.jl", 5)
convergence_test("elixir_euler_vortex_WF.jl", 5)
convergence_test("elixir_euler_vortex_FD.jl", 5)
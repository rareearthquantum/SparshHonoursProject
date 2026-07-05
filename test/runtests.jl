import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
using Test

# Include source files with paths anchored to `@__DIR__` as tests are added.

# Run test suite
println("Starting tests")

include("../src/input_pulse_methods.jl")
include("../src/echo_propagation.jl")

@testset "echo propagation stability" begin
    cfg = EchoConfig(Nd=8, Nt=128, Nz=16, alpha=100.0, d_width=200.0)
    result = run_propagation(cfg)

    @test result.atom_substeps == 8
    @test all(isfinite, result.Omega)
    @test all(isfinite, result.P)
end
ti = time()

@testset "SparshHonoursProject tests" begin
    @test 1 == 1
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")

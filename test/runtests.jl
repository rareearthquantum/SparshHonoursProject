import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
using Test
using Unitful

# Include source files with paths anchored to `@__DIR__` as tests are added.

# Run test suite
println("Starting tests")

include("../src/echo_propagation.jl")

@testset "automatic detuning grid" begin
    cfg = EchoConfig()
    dt = (cfg.Tf - cfg.Ti) / (cfg.Nt - 1)

    @test cfg.d_width ≈ 1 / dt
    @test cfg.Nd == ceil(Int, ND_TO_NT_RATIO * cfg.Nt)

    refined = EchoConfig(Nt=1000)
    @test refined.d_width ≈ (refined.Nt - 1) / (refined.Tf - refined.Ti)
    @test refined.Nd == 200
end

@testset "pulse-driven automatic config" begin
    pulses = PulseParams[
        PulseParams(center=1.0, width=0.2, area=pi/8),
        PulseParams(center=4.0, width=0.1, area=pi),
    ]
    cfg = config_from_pulses(pulses; samples_per_width=10)
    scaled = config_from_pulses([
        PulseParams(center=1e-6*p.center, width=1e-6*p.width, area=p.area)
        for p in pulses
    ]; samples_per_width=10)

    @test cfg.pulses == pulses
    @test cfg.Tf >= 2pulses[end].center - pulses[1].center
    @test cfg.Nt == scaled.Nt
    @test cfg.Nd == scaled.Nd
    @test scaled.d_width ≈ 1e6*cfg.d_width
    @test make_omega_input(cfg)(pulses[1].center) > 0

    physical = PulseParams(center=1u"μs", width=200u"ns", area=pi/8)
    @test physical.center == 1e-6
    @test physical.width == 200e-9
    @test config_from_pulses(pulses; Zf=2u"mm").Zf == 2e-3
    @test config_from_pulses(pulses; alpha=3u"s^-1*m^-1").alpha == 3.0
end

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

using Test
using Dates

include(joinpath(@__DIR__, "..", "src", "input_pulse_methods.jl"))
include(joinpath(@__DIR__, "..", "src", "config_hardcoded.jl"))
include(joinpath(@__DIR__, "..", "src", "integrators.jl"))
include(joinpath(@__DIR__, "..", "src", "simulation.jl"))
include(joinpath(@__DIR__, "..", "src", "saving.jl"))

@testset "Time coordinates in results and saved data" begin
    cfg = EchoConfig(Nt=32, Ti=2.0, Tf=4.0, Nd=1, Nz=1, Ny=1)
    result = run_2d_propagation(cfg)
    @test result.t_vec == make_time_grid(cfg)
    @test length(result.t_vec) == size(result.Omega, 1) == size(result.P, 1)
    mktempdir() do directory
        path = joinpath(directory, "result.jld2")
        save_result(path, result, 0.0)
        saved = JLD2.load(path)
        @test saved["t_vec"] == collect(result.t_vec)
        @test saved["Omega"] == result.Omega
        @test saved["P"] == result.P
    end
end

@testset "Pulse presets respect transverse width" begin
    for preset in (default_echo_2d_pulses, default_soliton_2d_pulses,
        default_4pi_2d_pulses, default_6pi_2d_pulses,
        default_weak_2d_pulses, default_pi_2d_pulses)
        pulses = preset(; y_pulse_width=0.3)
        @test all(pair -> pair[2].width == 0.3, pulses)
        center = pulses[1][1].center
        @test pulse_2d_sum(center, 0.3, pulses) / pulse_2d_sum(center, 0.0, pulses) ≈ exp(-0.5)
    end
end

@testset "Zero and one select the same single-point grids" begin
    reference = run_2d_propagation(EchoConfig(Nt=32, Nd=1, Nz=1, Ny=1))
    for Nd in (0, 1), Nz in (0, 1), Ny in (0, 1)
        cfg = EchoConfig(; Nt=32, Nd, Nz, Ny, d_width=0.0,
            Zi=2.0, Zf=2.0, Yi=3.0, Yf=3.0)
        @test (cfg.Nd, cfg.Nz, cfg.Ny) == (1, 1, 1)
        @test make_detunings(cfg) == [0.0]
        @test make_z_grid(cfg) == [0.0]
        @test make_y_grid(cfg) == [0.0]
        @test make_ky_grid(cfg) == [0.0]
        result = run_2d_propagation(cfg)
        @test size(result.Omega) == (32, 1, 1)
        @test size(result.P) == (32, 1)
        @test result.Omega == reference.Omega
        @test result.P == reference.P
    end

    cfg = EchoConfig(Nt=32, Nd=0, Nz=3, Ny=0, alpha=20.0)
    zero_result = run_2d_propagation(cfg)
    one_result = run_2d_propagation(EchoConfig(Nt=32, Nd=1, Nz=3, Ny=1, alpha=20.0))
    @test zero_result.Omega == one_result.Omega
    @test zero_result.P == one_result.P
    @test size(zero_result.Omega) == (32, 3, 1)
    @test all(isfinite, zero_result.P)

    cfg_2d = EchoConfig(Nt=32, Nd=0, Nz=3, Ny=4, alpha=0.0, beta=0.0)
    result_2d = run_2d_propagation(cfg_2d)
    @test size(result_2d.Omega) == (32, 3, 4)
    @test result_2d.Omega[:, end, :] ≈ result_2d.Omega[:, 1, :]
end

@testset "Invalid configurations fail early" begin
    for Nt in (-1, 0, 1)
        @test_throws ArgumentError EchoConfig(; Nt)
    end
    for name in (:Nd, :Nz, :Ny)
        @test_throws ArgumentError EchoConfig(; name => -1)
    end
    for Nz in (0, 1)
        @test_throws ArgumentError EchoConfig(; Nz, Ny=2)
    end
    @test_throws ArgumentError EchoConfig(Nd=2, d_width=0.0)
    @test_throws ArgumentError EchoConfig(Nd=0, d_width=-1.0)
    @test_throws ArgumentError EchoConfig(Nz=2, Zi=1.0, Zf=1.0)
    @test_throws ArgumentError EchoConfig(Ny=2, Yi=1.0, Yf=1.0)
    @test_throws ArgumentError EchoConfig(Ti=1.0, Tf=1.0)
    @test_throws ArgumentError EchoConfig(y_pulse_width=-1.0)
    for name in (:Ti, :Tf, :d_width, :Zi, :Zf, :alpha, :beta, :y_pulse_width, :Yi, :Yf)
        for value in (Inf, NaN)
            @test_throws ArgumentError EchoConfig(; name => value)
        end
    end
    valid = PulseParams(0.0, 1.0, 1.0)
    for invalid in (PulseParams(0.5, 0.0, pi), PulseParams(0.5, -0.1, pi),
        PulseParams(Inf, 0.1, pi), PulseParams(0.5, Inf, pi), PulseParams(0.5, 0.1, NaN))
        @test_throws ArgumentError EchoConfig(pulses=[(invalid, valid)])
        @test_throws ArgumentError EchoConfig(pulses=[(valid, invalid)])
    end
end

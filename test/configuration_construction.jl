# This test file can also run alone, without loading the solver or dependencies:
# julia --project=. test/configuration_construction.jl
using Test
include(joinpath(@__DIR__, "..", "src", "config", "simulation_config.jl"))

@testset "User-input configuration construction" begin
    @testset "Defaults and inspection" begin
        cfg = SimulationConfig()
        @test cfg.geometry == :two_d
        @test cfg.pulses isa Vector{PulseSpec}
        @test length(cfg.pulses) == 1
        @test cfg.pulses[1].shape == :gaussian
        @test cfg.pulses[1].center == 0.5
        @test cfg.pulses[1].width == 0.1
        @test cfg.pulses[1].area ≈ pi/8
        @test cfg.pulses[1].y_center == 0.0
        @test cfg.pulses[1].y_width == 1.0
        @test (cfg.medium.d_width, cfg.medium.alpha, cfg.medium.beta) == (1000.0, 0.0, 1.0)
        @test (cfg.grid.Nt, cfg.grid.Nd, cfg.grid.Nz, cfg.grid.Ny) == (256, 128, 64, 128)
        @test (cfg.grid.Ti, cfg.grid.Tf, cfg.grid.Zi, cfg.grid.Zf) == (0.0, 1.0, 0.0, 1.0)
        @test (cfg.grid.Yi, cfg.grid.Yf) == (-3.0, 3.0)
        for value in (cfg, cfg.pulses[1], cfg.medium, cfg.grid)
            @test occursin(string(nameof(typeof(value))), sprint(show, value))
        end
    end

    @testset "Explicit physical inputs" begin
        # Example convention: seconds, metres, angular frequencies in rad/s.
        pulse = PulseSpec(shape=:top_hat, center=2e-6, width=0.2e-6,
            area=pi, y_center=0.1e-3, y_width=0.5e-3)
        medium = MediumSpec(d_width=2pi*1e6, alpha=2e8, beta=1e-6)
        grid = GridSpec(Nt=64, Ti=1e-6, Tf=4e-6, Nd=16, Nz=8,
            Zi=0, Zf=0.01, Ny=1, Yi=-0.002, Yf=0.002)
        cfg = SimulationConfig(geometry=:one_d, pulses=[pulse], medium=medium, grid=grid)
        @test cfg.pulses[1] === pulse
        @test cfg.medium === medium
        @test cfg.grid === grid
        @test cfg.geometry == :one_d
        @test (pulse.center, pulse.width, pulse.area) == (2e-6, 0.2e-6, Float64(pi))
        @test (pulse.shape, pulse.y_center, pulse.y_width) == (:top_hat, 0.1e-3, 0.5e-3)
        @test (medium.d_width, medium.alpha, medium.beta) == (2pi*1e6, 2e8, 1e-6)
        @test (grid.Nt, grid.Nd, grid.Nz, grid.Ny) == (64, 16, 8, 1)
        @test (grid.Ti, grid.Tf, grid.Zi, grid.Zf, grid.Yi, grid.Yf) ==
            (1e-6, 4e-6, 0.0, 0.01, -0.002, 0.002)
    end

    @testset "No implicit resolution or shared defaults" begin
        # Labels and window choices are stored only; they do not expand recipes
        # or silently rewrite explicit counts, pulse widths, or transverse bounds.
        grid = GridSpec(Ti=2.0, Tf=4.0, Nd=0, Nz=0, Ny=0, Yi=7.0, Yf=9.0)
        cfg = SimulationConfig(geometry=:zero_d, grid=grid)
        @test (cfg.grid.Nd, cfg.grid.Nz, cfg.grid.Ny) == (0, 0, 0)
        @test (cfg.grid.Yi, cfg.grid.Yf) == (7.0, 9.0)
        @test (cfg.pulses[1].center, cfg.pulses[1].width) == (0.5, 0.1)
        @test isempty(SimulationConfig(pulses=PulseSpec[]).pulses)
        sequence = [PulseSpec(area=pi/2), PulseSpec(center=0.8, area=pi)]
        @test SimulationConfig(pulses=sequence).pulses == sequence
        first_cfg, second_cfg = SimulationConfig(), SimulationConfig()
        @test first_cfg.pulses !== second_cfg.pulses
        push!(first_cfg.pulses, PulseSpec())
        @test length(second_cfg.pulses) == 1
    end
end

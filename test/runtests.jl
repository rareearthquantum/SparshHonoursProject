using DrWatson, Test
@quickactivate "SparshHonoursProject"

include(srcdir("simulation_config.jl"))

# Run test suite
println("Starting tests")
ti = time()

@testset "SparshHonoursProject tests" begin
    dimensionless = load_simulation_config(
        scriptsdir("input_params_large_single_pulse.toml")
    )
    @test dimensionless.mode == :dimensionless
    @test dimensionless.N == (128, 64, 64, 64)
    @test dimensionless.alpha == 0.1
    @test dimensionless.beta == 0.02
    @test length(dimensionless.pulse_params) == 1

    physical = load_simulation_config(
        scriptsdir("input_params_2_10micron_y_pulses_40micron_seperation.toml")
    )
    @test physical.mode == :physical_si
    @test physical.alpha ≈ 0.1
    @test physical.beta ≈ 2.3856725796184706
    @test physical.ranges == ((0.0, 1.0), (0.0, 5.0), (-5.0, 5.0))
    @test length(physical.pulse_params) == 2
    @test physical.pulse_params[1][2].center == -2.0
    @test physical.pulse_params[2][2].center == 2.0
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")

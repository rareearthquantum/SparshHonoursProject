# Fast physics regressions, loaded by runtests.jl after the solver components.
# All values use a consistent numerical unit system; no plots or files are made.

@testset "Local constant-drive Rabi evolution" begin
    for Nd in (1, 5)
        cfg = EchoConfig(Nt=129, Ti=0.3, Tf=1.3, Nd=Nd, d_width=4.0,
            Nz=1, Ny=1, alpha=0.0, beta=0.0)
        drive = Float64(pi)
        result = run_2d_propagation(cfg;
            omega_2d_input=(t, y) -> drive, show_progress=false)
        elapsed = result.t_vec .- cfg.Ti

        # Exact coherence for a constant real drive, ground-state initial atoms,
        # and angular detuning delta. It is in the unrotated frame returned as P.
        # For delta=0 this reduces to (im/2)*sin(drive*elapsed), a pi rotation.
        expected = map(elapsed) do tau
            sum(result.detunings) do delta
                frequency = hypot(drive, delta)
                drive * delta / (2frequency^2) * (1 - cos(frequency * tau)) +
                    im * drive / (2frequency) * sin(frequency * tau)
            end / Nd
        end

        @test size(result.Omega) == (129, 1, 1)
        @test result.Omega[:, 1, 1] == fill(drive, cfg.Nt)
        # dt=1/128 resolves the largest generalised Rabi angle per step (<0.03).
        # 2e-8 allows accumulated RK4 truncation, well below the coherence scale 1/2.
        @test maximum(abs.(result.P[:, 1] .- expected)) <= 2e-8
        @test maximum(abs, result.P) <= 0.5 + 2e-8
    end

    cfg = EchoConfig(Nt=17, Nd=5, d_width=4.0, Nz=1, Ny=1)
    dark = run_2d_propagation(cfg;
        omega_2d_input=(t, y) -> 0.0, show_progress=false)
    @test all(iszero, dark.Omega)
    @test all(iszero, dark.P)
end

@testset "1D weak-field propagation" begin
    drive = 1e-3
    errors = Float64[]
    for Nz in (9, 17)
        cfg = EchoConfig(Nt=65, Ti=0.0, Tf=1.0, Nd=1, d_width=0.0,
            Nz=Nz, Zi=0.0, Zf=0.5, Ny=1, alpha=0.4, beta=0.0)
        result = run_2d_propagation(cfg;
            omega_2d_input=(t, y) -> drive, show_progress=false)

        # At resonance and weak drive, inversion remains approximately -1:
        # dOmega/dz = -(alpha/2)*integral_0^t Omega(s,z) ds.
        # A constant boundary drive has solution drive*J0(sqrt(2alpha*z*t)).
        # Evaluate its power series directly, avoiding a special-function dependency.
        # x<=0.1: the omitted term after k=6 is <4e-15 relative to the drive.
        expected = [drive * sum((-cfg.alpha * z * t / 2)^k / factorial(k)^2
            for k in 0:6) for t in result.t_vec, z in result.z_vec]
        field = result.Omega[:, :, 1]
        error = maximum(abs.(field .- expected))
        push!(errors, error)

        @test size(result.Omega) == (65, Nz, 1)
        @test all(isfinite, result.P)
        @test field[:, 1] == fill(drive, cfg.Nt)
        @test maximum(abs, imag.(field)) <= 1e-14
        @test real(field[end, end]) < drive
        # Euler startup + AB2 gives O(dz^2) global error for this smooth solution.
        # Bounds allow that error and the O(drive^3) nonlinear correction; these
        # coarse grids are regression checks, not scientific convergence claims.
        @test error <= (Nz == 9 ? 6e-8 : 2e-8)
    end
    @test errors[2] < errors[1] / 3
end

@testset "2D exact diffraction and transverse norm" begin
    Ny = 8
    # Eight distinct samples of a 2pi-periodic domain, with no duplicate endpoint.
    # The inclusive solver grid then has spacing 2pi/Ny and FFT modes k=1,2.
    cfg = EchoConfig(Nt=17, Nd=3, d_width=2.0, Nz=5, Zi=0.2, Zf=0.8,
        Ny=Ny, Yi=-pi, Yf=pi - 2pi/Ny, alpha=0.0, beta=0.3)
    amplitude(t) = 0.25 + 0.1t
    input(t, y) = amplitude(t) * (1 + 0.3cos(y) + 0.2sin(2y))
    result = run_2d_propagation(cfg; omega_2d_input=input, show_progress=false)

    # With alpha=0, each Fourier mode gains exp(-im*beta*k^2*(z-Zi)).
    # The superposition checks diffraction phase, sign, and k^2 dependence.
    expected = [amplitude(t) * (1 +
        0.3cos(y) * cis(-cfg.beta * (z-cfg.Zi)) +
        0.2sin(2y) * cis(-4cfg.beta * (z-cfg.Zi)))
        for t in result.t_vec, z in result.z_vec, y in result.y_vec]
    @test size(result.Omega) == (17, 5, 8)
    @test all(isfinite, result.P)
    # Diffraction is an exact spectral phase update here; only roundoff accumulates.
    @test maximum(abs.(result.Omega .- expected)) <= 3e-14
    @test maximum(abs.(result.Omega[:, end, :] .- result.Omega[:, 1, :])) > 1e-2
    energy = dropdims(sum(abs2.(result.Omega); dims=3); dims=3)
    @test maximum(abs.(energy .- energy[:, 1])) <= 1e-13
end

@testset "Uniform transverse field reduces coupled 2D to 1D" begin
    common = (; Nt=65, Nd=3, d_width=2.0, Nz=5, Zi=0.0, Zf=0.25,
        alpha=0.4, beta=0.3)
    input(t, y) = 0.2 + 0.03t
    one_d = run_2d_propagation(EchoConfig(; common..., Ny=1);
        omega_2d_input=input, show_progress=false)
    two_d = run_2d_propagation(EchoConfig(; common..., Ny=8, Yi=-pi, Yf=3pi/4);
        omega_2d_input=input, show_progress=false)

    # A uniform transverse input occupies only k=0: diffraction vanishes even
    # though beta is nonzero, while atom-field coupling is active in both runs.
    @test maximum(abs.(one_d.Omega[:, end, :] .- one_d.Omega[:, 1, :])) > 1e-4
    for y in axes(two_d.Omega, 3)
        # Identical discrete equations; allow accumulated FFT/Float64 roundoff.
        @test maximum(abs.(two_d.Omega[:, :, y] .- one_d.Omega[:, :, 1])) <= 1e-12
        @test maximum(abs.(two_d.P[:, y] .- one_d.P[:, 1])) <= 1e-12
    end
end

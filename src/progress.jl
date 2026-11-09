"""
    SimulationProgress(total_steps; enabled=true, io=stdout)

Report actual completion at roughly ten-percent intervals. A work step is one
propagation step or the final polarisation pass. Elapsed time uses a monotonic
clock, and each report occupies a new line for readable batch-job logs.
"""
mutable struct SimulationProgress{I<:IO}
    total_steps::Int
    next_percent::Int
    started_ns::UInt64
    enabled::Bool
    io::I
end

function SimulationProgress(total_steps::Integer; enabled::Bool=true, io::IO=stdout)
    total_steps > 0 || throw(ArgumentError("total_steps must be positive"))
    return SimulationProgress(Int(total_steps), 10, time_ns(), enabled, io)
end

function report_progress!(progress::SimulationProgress, completed_steps::Integer)
    progress.enabled || return nothing
    0 <= completed_steps <= progress.total_steps ||
        throw(ArgumentError("completed_steps must be between 0 and total_steps"))

    percent = fld(100 * completed_steps, progress.total_steps)
    percent >= progress.next_percent || return nothing

    elapsed = round((time_ns() - progress.started_ns) / 1e9; digits=2)
    println(progress.io,
        "$(percent)% complete ($(completed_steps)/$(progress.total_steps) steps) in $(elapsed)s")
    flush(progress.io)
    progress.next_percent = 10 * (percent ÷ 10 + 1)
    return nothing
end

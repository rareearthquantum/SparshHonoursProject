function s_evolve!(du, u, p, t)

    d, a_interp = p

    du .= im .* d .* u + im .* a_interp(t) .* ones(length(d))
end

function f_a_2d!(alpha::Float64, s_grid::AbstractArray, beta::Float64, dy2::ComplexF64)
    im * (alpha * sum(s_grid) + beta * dy2)
end


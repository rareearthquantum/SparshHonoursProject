function s_evolve!(du, u, p, t)

    d, a_interp = p

    du .= im .* d .* u .+ im .* a_interp(t) 
end

function f_a_2d!(alpha::Real, s_grid::AbstractArray, beta::Real, dy2::Complex)
    im * (alpha * sum(s_grid) + beta * dy2)
end


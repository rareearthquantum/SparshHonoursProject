include(srcdir("stepping_methods.jl"))

function stepping_ab_4step!(v::AbstractArray, f::Function, h::Float64, i::Int64)
    (i >= length(v)) && return
    if (i >= 4)
        adams_bashforth_4step!(v, f, h, i)
    else
        if (i == 1)
            euler_step!(v, f, h, i)
        elseif (i == 2)
            adams_bashforth_2step!(v, f, h, i)
        elseif (i == 3)
            adams_bashforth_3step!(v, f, h, i)
        end
    end
end

function stepping_ab_3step!(v::AbstractArray, f::Function, h::Float64, i::Int64)
    (i >= length(v)) && return
    if (i >= 3)
        adams_bashforth_3step!(v, f, h, i)
    else
        if (i == 1)
            euler_step!(v, f, h, i)
        elseif (i == 2)
            adams_bashforth_2step!(v, f, h, i)
        end
    end
end

function stepping_ab_2step!(v::AbstractArray, f::Function, h::Float64, i::Int64)
    (i >= length(v)) && return
    if (i >= 2)
        adams_bashforth_2step!(v, f, h, i)
    else
        euler_step!(v, f, h, i)
    end
end

include("../constants.jl")

const frequency = 200tera
const omega0 = 2.0 * pi * frequency
const c = 3e8
const k = omega0 / c
const mu0 = 1.25663706e-6

const beta = 1/(2k)
@show beta
const alpha = 0.0 #mu0*omega0/2
@show alpha

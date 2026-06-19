const tera, giga, mega, kilo, centi, milli, micro, nano, pico = 1e12, 1e9, 1e6, 1e3, 1e-2, 1e-3, 1e-6, 1e-9, 1e-12

const frequency = 200tera
const omega0 = 2.0 * pi * frequency
const c = 3e8
const k = omega0 / c

const beta = 1/(2k) * 1/(2pi)

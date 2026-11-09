# Parameter labels shared by data and figure filenames.
"""
    make_param_info(cfg)

Return a filename label containing grid sizes, detuning width, coupling,
diffraction, pulse count, and each temporal pulse area and width.
"""
function make_param_info(cfg)
    parameter_info = "_Nt=$(cfg.Nt)_Nd=$(cfg.Nd)_Nz=$(cfg.Nz)_Ny=$(cfg.Ny)_dwidth=$(cfg.d_width)_alpha=$(cfg.alpha)_beta=$(cfg.beta)_pulsecount=$(length(cfg.pulses))_"
    for pulse in cfg.pulses
        parameter_info *= "_area=$(pulse[1].area)_width=$(pulse[1].width)_"
    end
    return parameter_info
end

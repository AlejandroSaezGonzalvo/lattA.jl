module lattA

import juobs, ADerrors

include("types.jl")
export EnsInfo

include("reader.jl")
export get_corr_wil, get_corr_tm, get_dSdm, get_YW

include("obs.jl")
export get_m, get_mpcac, get_f_wil, get_f_tm, get_t0

include("tools.jl")
export model_av, pvalue, get_model, fve, fit_alg, fit_defs

end

module lattA

import juobs, ADerrors

include("types.jl")
export EnsInfo

include("reader.jl")
export get_corr_wil, get_corr_tm, get_dSdm, get_YW, read_ens_wil, read_ens_tm, read_ens_tm_sym, read_ens_csv

include("obs.jl")
export get_m, get_mpcac, get_f_wil, get_f_tm, get_t0

include("tools.jl")
export model_av, pvalue, get_model, fve, fit_alg, fit_defs

include("fit_funs.jl")
export md_1, md_2, match_sym, match_constTr, interp_fpik_sym, interp_fpik_constTr

include("plot.jl")
export matching_sym_plot, interp_fpik_sym_plot

include("csv_tools.jl")
export csv2Corr

end

module lattA

using ADerrors: err

import juobs, ADerrors

include("types.jl")
export EnsInfo

include("reader.jl")
export get_corr_wil, get_corr_tm, get_dSdm, get_YM, read_ens_wil, read_ens_tm, read_ens_tm_sym, read_ens_csv, concat_data, read_mesons_correction_multichunks, read_mesons_multichunks, get_corr_TSM_multichunks, get_corr_TSM, read_ens_TSM

include("obs.jl")
export get_m, get_m_pbc, get_mpcac, get_f_wil, get_f_wil_pbc, get_f_tm, get_f_tm_pbc, get_t0, get_w0t0
export get_m_ALPHA, get_m_pbc_ALPHA, get_mpcac_ALPHA, get_f_wil_ALPHA, get_f_tm_ALPHA, get_t0_ALPHA

include("tools.jl")
export model_av, pvalue, get_model, fve, fve_inv, fit_alg, fit_alg_LBFGS, fit_defs, corr_sym_E250, corr_sym_D450

include("fit_funs.jl")
export md_1, md_2, match_m12_sym, match_phi2_sym, match_m12, match_phi2, match_phi4, interp_fpik_sym, interp_fpik_constTr

include("csv_tools.jl")
export csv2Corr

end

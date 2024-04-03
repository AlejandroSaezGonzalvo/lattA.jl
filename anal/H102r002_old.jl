A = 3
using Revise, juobs, BDIO, ADerrors, PyPlot, Gnuplot, LsqFit, Distributions, Dates, SpecialFunctions

cd("/home/asaez/cls_ens/codes/analysis_cls/PDG")
include("/home/asaez/cls_ens/codes/analysis_cls/csv_tools.jl") ; include("/home/asaez/cls_ens/codes/analysis_cls/analysis_st_syst_tmin_tmax_1q_combined.jl") ; include("/home/asaez/cls_ens/codes/analysis_cls/analysis_tm_syst_tmin_tmax_1q_combined.jl") ; include("/home/asaez/cls_ens/codes/analysis_cls/latex_table.jl") ; include("/home/asaez/cls_ens/codes/analysis_cls/paths_combined.jl") ; include("/home/asaez/cls_ens/codes/analysis_cls/PDG/inputs_1098.jl") ; include("/home/asaez/cls_ens/codes/analysis_cls/fit_models.jl") ; include("/home/asaez/cls_ens/codes/analysis_cls/fit_routine_corr.jl") ; include("/home/asaez/cls_ens/codes/analysis_cls/plateaus.jl") 

#include("/home/asaez/cls_ens/codes/analysis_cls/obs_ethan+algebra_syst_tmin_tmax_array.jl")
include("/home/asaez/cls_ens/codes/analysis_cls/juobs_tools_new.jl")
cd("/home/asaez/cls_ens/codes/analysis_cls")

d = Normal(0.0,10.0)

global Vandermonde = 1 ## if Vandermonde = 1, uses Vandermonde & Bayes, otherwise plateaus by hand and no Vandermonde
A = parse(Int64,ARGS[1])
#A = 1
chosen_ens = A:A 

global t0_bayesian = true ## if false, computes t0 with hand plateau
improved = 1
if improved == 0 ## 1 if mpcac improved with cA, 0 if unimproved with cA
	ca .= 0
end

#================================ choose and analyze one ensemble from complete list =====================================================#

	global ens = ["H101", "H102r001", "H102r002", "H105", "H105r005", "H400", "N202", "N203", "N200", "D200", "N300r002", "N302", "J303"]
	global ens_tm = ens
	global sym_ens = ["H101", "H400", "N202", "N300r002"]

	## choose the observables you are interested in:
		ens = ens[chosen_ens]
		ens_tm = ens_tm[chosen_ens]
		path = path[chosen_ens]
		path_tm = path_tm[chosen_ens]
		path_ll = path_ll[chosen_ens]
		path2_ll = path2_ll[chosen_ens]
		path_ss = path_ss[chosen_ens]
		path2_ss = path2_ss[chosen_ens]
		path_ls = path_ls[chosen_ens]
		path2_ls = path2_ls[chosen_ens]
		path_ll_sym = path_ll_sym[chosen_ens]
		path2_ll_sym = path2_ll_sym[chosen_ens]
		path2_ls_sym = path2_ls_sym[chosen_ens]
		path_ss_sym = path_ss_sym[chosen_ens]
		path2_ss_sym = path2_ss_sym[chosen_ens]
		path_ls_sym = path_ls_sym[chosen_ens]
		path_ms = path_ms[chosen_ens]
		sdr_path = sdr_path[chosen_ens]
		path_rwf = path_rwf[chosen_ens]
		L = L[chosen_ens]
		T = T[chosen_ens]
		beta = beta[chosen_ens]
		delta_b = delta_b[chosen_ens]
		bAtil = bAtil[chosen_ens]
		bPtil = bPtil[chosen_ens]
		ZA = ZA[chosen_ens]
		ca = ca[chosen_ens]
		mul_guess = mul_guess[chosen_ens]
		mus_guess = mus_guess[chosen_ens]
		plat_t0 = plat_t0[chosen_ens] ; plat_mpi = plat_mpi[chosen_ens] ; plat_mk = plat_mk[chosen_ens] ; plat_fpi = plat_fpi[chosen_ens] ; plat_fk = plat_fk[chosen_ens] ; plat_12 = plat_12[chosen_ens] ; plat_13 = plat_13[chosen_ens] ; plat_mpi_tm = plat_mpi_tm[chosen_ens] ; plat_mk_tm = plat_mk_tm[chosen_ens] ; plat_fpi_tm = plat_fpi_tm[chosen_ens] ; plat_fk_tm = plat_fk_tm[chosen_ens] ; plat_12_tm = plat_12_tm[chosen_ens]
		#(mpi,mk,m12,m13,Rpi,Rk)_guess
	##

	matching_plots_sym = ["H101"] 
	matching_plots = ["H105"] 

	wpm = Dict{String, Vector{Float64}}()
	for INDEX in 1:length(ens)
		global wpm[ens[INDEX]] = wpm[ens_tm[INDEX]] = [-1.0, -1.0, 4.0, -1.0]
	end

	start_time = Dates.format(now(), "yyyy-mm-dd at HH:MM")
	println("Start: ", start_time)



using juobs, ADerrors, CSV, DataFrames

@doc raw"""
    csv2Corr(path::String, ivrep::Union{Nothing, Vector{Int64}}=nothing; trunc::Union{Nothing, Vector{Int64}}=nothing, id::Union{Nothing, String}=nothing)     

Reads a csv file produced by xml2csv.m (MATLAB) and returns a juobs.Corr. Values of Gamma, kappa and mu are determined using the filename.
The vector `ivrep` has to be fixed to number of configurations of each replica for files 
that contain two or more replicas. The extra argument `trunc` allows to truncate each replica. `id` can be specified as an extra argument.
```@example
corr = csv2Corr(path_to_csv)
corr = csv2Corr(path_to_csv, [505, 540])
corr_trunc = csv2Corr(path_to_csv, [505, 540], trunc=[450, 540])
corr_trunc = csv2Corr(path_to_csv, [505, 540], trunc=[450, 540], id="H400")
```
"""

function csv2Corr(path::String, ivrep::Union{Nothing, Vector{Int64}}=nothing; trunc::Union{Nothing, Vector{Int64}}=nothing, id::Union{Nothing, String}=nothing, flag_strange::Bool=false)
    gamma_dict = Dict([("PS", "G5"), ("A0", "G0G5"), ("G0", "G0")])
    bname = basename(path)
    
    #REGEX y0
    x = match(r"s[0-9]+", bname)
    y0 = parse(Int64, x.match[2:end])
    
    #REGEX kappa1, kappa2
    k = Vector{Float64}(undef, 0)
    for x in eachmatch(r"k[0-9]\.[0-9]+", bname)
        aux = parse(Float64, x.match[2:end])
        push!(k, aux)
    end
    #REGEX mu1, mu2
    mu = Vector{Float64}(undef, 0)
    for x in eachmatch(r"mup[0-9]\.[0-9]+", bname)
        aux = parse(Float64, x.match[4:end])
        push!(mu, aux)
    end
    #REGEX ensemble id
    if isnothing(id)
        x = match(r"[A-Z][0-9]{3}r?[0-9]?[0-9]?[0-9]?", bname)
        id = convert(String, x.match)
    end
    #REGEX Gamma1, Gamma2
    gamma = Vector{String}(undef, 0)
    g1, g2 = bname[1:2], bname[3:4]
    push!(gamma, gamma_dict[g1])
    push!(gamma, gamma_dict[g2])

    #READ csv
    df  = CSV.read(path, DataFrame, header=false)
    data = Matrix(df)

    if id == "H105r005" && flag_strange == true
        flagged_cfg = [254, 255, 256, 257, 259, 260, 261, 264, 265, 266, 269, 280, 282, 283, 284, 285, 286, 287, 288, 289, 291, 299, 301, 313, 314, 315, 316, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342]
        W_H105r005 = [1.0 for i in 1:length(data[:,1])]
        for i in flagged_cfg
            W_H105r005[i] = -1.0
        end
        data = data .* W_H105r005
        W_obs = uwreal(W_H105r005, id)
    end
    
    #Checks + truncation
    if isnothing(ivrep)
        ivrep = [size(data, 1)]
    end
    if sum(ivrep) != size(data, 1)
        error("Dimension mismatch: ivrep")
    end

    if isnothing(trunc)
        aux_data = data
        aux_ivrep = ivrep
    else
        if any(trunc .> ivrep)
            error("Dimension mismatch: trunc")
        end
        aux_data = data[1:trunc[1], :]
        ptr_cnfg = ivrep[1]
        for k = 2:length(trunc)
            aux_data = vcat(aux_data, data[ptr_cnfg + 1:ptr_cnfg + trunc[k], :])
            ptr_cnfg = ptr_cnfg + ivrep[k]
        end
        aux_ivrep = trunc
    end

    if flag_strange == true && id == "H105r005" 
        ow = [uwreal(aux_data[:, k], id, aux_ivrep) for k = 1:size(data, 2)]
        obs = [ow[x0] / W_obs for x0 = 1:size(data, 2)]
    else
        obs = [uwreal(aux_data[:, k], id, aux_ivrep) for k = 1:size(data, 2)]
    end    

    return juobs.Corr(obs, k, mu, gamma, y0, [.0,.0,.0], [.0,.0,.0])
end

err = ADerrors.err

INDEX = 1

global inner_index = INDEX

		println("\nAnalyzing standard observables : ", ens[INDEX])
		
        (dm, dmpi_s1, dmpi_s2, dmk_s1, dmk_s2, dm12_s1, dm12_s2, dfpi_s1, dfpi_s2, dfk_s1, dfk_s2, dmpisyst_s1, dmpisyst_s2, dmksyst_s1, dmksyst_s2, dm12syst_s1, dm12syst_s2, dRpisyst_s1, dRpisyst_s2, dRksyst_s1, dRksyst_s2, dm13_s1, dm13_s2, dm34_s1, dm34_s2, dm13syst_s1, dm13syst_s2, dm34syst_s1, dm34syst_s2) = analyze_st_2(path[INDEX], path_rwf[INDEX], path_ms[INDEX], sdr_path[INDEX], ens[INDEX], INDEX) 
		derivatives = [dmpi_s1, dmpi_s2, dmk_s1, dmk_s2, dm12_s1, dm12_s2, dfpi_s1, dfpi_s2, dfk_s1, dfk_s2, dmpisyst_s1, dmpisyst_s2, dmksyst_s1, dmksyst_s2, dm12syst_s1, dm12syst_s2, dRpisyst_s1, dRpisyst_s2, dRksyst_s1, dRksyst_s2, dm13_s1, dm13_s2, dm34_s1, dm34_s2, dm13syst_s1, dm13syst_s2, dm34syst_s1, dm34syst_s2]
			
		
        ## standard analysis
		uwerr(t0_sh[INDEX]) ; uwerr(mpi_sh[INDEX]) ; uwerr(mk_sh[INDEX]) ; uwerr(m12_I_sh[INDEX]) ; uwerr(m13_I_sh[INDEX]) ; uwerr(fpi_sh[INDEX]) ; uwerr(fk_sh[INDEX]) 
		println("\nt0_sh = ", t0_sh[INDEX], "\nmpi_sh = ", mpi_sh[INDEX], "\nmk_sh = ", mk_sh[INDEX], "\nm12_I_sh = ", m12_I_sh[INDEX], "\nm13_I_sh = ", m13_I_sh[INDEX], "\nm34_I_sh = ", m34_I_sh[INDEX], "\nfpi_sh = ", fpi_sh[INDEX], "\nfk_sh = ", fk_sh[INDEX])
		GC.gc()

		global dm = dm

        index = INDEX = 1
        identity = ens[index]
        path_corr = path[index]
        path_tm = path_tm[index]
        path_ms = path_ms[index]
        path_rwf = path_rwf[index]
        sdr_path = sdr_path[index]
        path_ll = path_ll[index]
        path_ls = path_ls[index]
        path_ss = path_ss[index]
        path_ll_sym = path_ll_sym[index]
        path_ls_sym = path_ls_sym[index]
        path_ss_sym = path_ss_sym[index]
        path2_ll = path2_ll[index]
        path2_ss = path2_ss[index]
        path2_ll_sym = path2_ll_sym[index]
        path2_ls = path2_ls[index]
        path2_ls_sym = path2_ls_sym[index]
        path2_ss_sym = path2_ss_sym[index]

data11 = csv2Corr.(path_ll, id=identity, flag_strange=true) 
				data12 = csv2Corr.(path_ls, id=identity, flag_strange=true)
				data212 = csv2Corr.(path2_ls, id=identity, flag_strange=true)
				data22 = csv2Corr.(path2_ll, id=identity, flag_strange=true) 
				data11_sym = csv2Corr.(path_ll_sym, id=identity, flag_strange=true) 
				data12_sym = csv2Corr.(path_ls_sym, id=identity, flag_strange=true)
				data212_sym = csv2Corr.(path2_ls_sym, id=identity, flag_strange=true)
				data22_sym = csv2Corr.(path2_ll_sym, id=identity, flag_strange=true) 
				data_strange = csv2Corr.(path_ss, id=identity, flag_strange=true)
				data_strange_sym = csv2Corr.(path_ss_sym, id=identity, flag_strange=true)
				data2_strange = csv2Corr.(path2_ss, id=identity, flag_strange=true)
				data2_strange_sym = csv2Corr.(path2_ss_sym, id=identity, flag_strange=true)


				data11_S = [corr_sym(data11[i],data11_sym[i],+1) for i in 1:1:length(data11)]
				data12_S = [corr_sym(data12[i],data12_sym[i],+1) for i in 1:1:length(data12)]
				data22_S = [corr_sym(data22[i],data22_sym[i],-1) for i in 1:1:length(data22)]
				data212_S = [corr_sym(data212[i],data212_sym[i],-1) for i in 1:1:length(data212)]
				data_strange_S = [corr_sym(data_strange[i],data_strange_sym[i],+1) for i in 1:1:length(data_strange)]
				data2_strange_S = [corr_sym(data2_strange[i],data2_strange_sym[i],-1) for i in 1:1:length(data2_strange)]

				k_tm = [data11[i].kappa for i in 1:length(data11)]
				kl_tm = [k_tm[i][1] for i in 1:length(k_tm)] ## the kappas are always degenerate in the sea in tmQCD ???
				ks_tm = [k_tm[i][2] for i in 1:length(k_tm)]
				mu_tm = [data11[i].mu for i in 1:length(data11)]
				mul_tm = [mu_tm[i][1] for i in 1:length(mu_tm)]

				global x1 = [kl_tm mul_tm mul_tm]
				global x2 = [kl_tm mul_tm mul_tm]

				k_tm = [data12[i].kappa for i in 1:length(data12)]
				kl_tm = [k_tm[i][1] for i in 1:length(k_tm)] ## the kappas are always degenerate in the sea in tmQCD ???
				ks_tm = [k_tm[i][2] for i in 1:length(k_tm)]
				mu_tm = [data12[i].mu for i in 1:length(data12)]
				mul_tm = [mu_tm[i][1] for i in 1:length(mu_tm)]
				mus_tm = [mu_tm[i][2] for i in 1:length(mu_tm)]

				global x3 = [kl_tm mul_tm mus_tm] ## x = [kl_tm mul_tm mus_tm]
				global x = [x1; x2; x3]

				#===============================================================================================================================#

				phi4_tm = Array{uwreal,1}()
				mpi_tm = Array{uwreal,1}()
				syst_mpi = Array{uwreal,1}()
				mk_tm = Array{uwreal,1}()
				syst_mk = Array{uwreal,1}()
				fk_tm = Array{uwreal,1}()
				syst_Rk = Array{uwreal,1}()
				m12_tm = Array{uwreal,1}()
				m13_tm = Array{uwreal,1}()
				m34_tm = Array{uwreal,1}()
				syst_m12 = Array{uwreal,1}()
				syst_m13 = Array{uwreal,1}()
				syst_m34 = Array{uwreal,1}()
				fpi_tm = Array{uwreal,1}()
				syst_Rpi = Array{uwreal,1}()
				mu1_pi = Array{Float64,1}()
				mu2_pi = Array{Float64,1}()
				mu1_k = Array{Float64,1}()
				mu2_k = Array{Float64,1}()						

				for i in 1:9
					println(i)
					(aux1, aux_syst) = m_ground(data11_S[i];plot_data=false,label=string("meff-pion-",i,"-",identity), wpm=wpm)
					(aux, aux_syst_12) = mpcac_bayesian(data22_S[i],data11_S[i];ca=ca[index],plot_data=false,label=string("mpcac-12-",i,"-",identity), wpm=wpm)
					push!(m12_tm, aux) ; push!(syst_m12, aux_syst_12)
					push!(mpi_tm, aux1) ; push!(syst_mpi, aux_syst)
					(aux, aux_syst) = dec_const_pcvc_bayesian(data11[i],data11_sym[i],aux1;plot_data=false,label=string("R-pion-",i,"-",identity), wpm=wpm)
					push!(fpi_tm, aux) ; push!(syst_Rpi, aux_syst)
					push!(mu1_pi, data11[i].mu[1])
					push!(mu2_pi, data11[i].mu[2])
					(aux, aux_syst_34) = mpcac_bayesian(data2_strange_S[i],data_strange_S[i];ca=ca[index],plot_data=false,label=string("mpcac-34-",i,"-",identity), wpm=wpm)
					push!(m34_tm, aux) ; push!(syst_m34, aux_syst_34)
				end

				for i in 1:27
					(aux2, aux_syst) = m_ground(data12_S[i];plot_data=false,label=string("meff-kaon-",i,"-",identity), wpm=wpm)
					uwerr(aux2)
					push!(mk_tm, aux2)
					push!(syst_mk, aux_syst)
					(aux, aux_syst) = dec_const_pcvc_bayesian(data12[i],data12_sym[i],aux2;plot_data=false,label=string("R-kaon-",i,"-",identity), wpm=wpm)
					push!(fk_tm, aux) ; push!(syst_Rk, aux_syst)
					push!(mu1_k, data12[i].mu[1])
					push!(mu2_k, data12[i].mu[2])
					(aux, aux_syst_13) = mpcac_bayesian(data212_S[i],data12_S[i];ca=ca[index],plot_data=false,label=string("mpcac-13-",i,"-",identity), wpm=wpm)
					push!(m13_tm, aux) ; push!(syst_m13, aux_syst_13)
				end

				phi2_tm = [8*t0[index] *mpi_tm[i] ^2 for i in 1:length(mpi_tm)]
				
				for i in 1:length(x1[:,1])
					for j in 1:length(x2[:,1])
						if x3[j,2] == x1[i,2]
							push!(phi4_tm, 8*t0[index] *(0.5 *mpi_tm[i] ^2 + mk_tm[j] ^2))
						end
					end
				end

                for i in 1:length(mpi_tm) uwerr(mpi_tm[i],wpm) end
			for i in 1:length(mk_tm) uwerr(mk_tm[i],wpm) end
			for i in 1:length(m12_tm) uwerr(m12_tm[i],wpm) end
			for i in 1:length(m13_tm) uwerr(m13_tm[i],wpm) end
			for i in 1:length(m34_tm) uwerr(m34_tm[i],wpm) end
			for i in 1:length(phi2_tm) uwerr(phi2_tm[i],wpm) end
			for i in 1:length(phi4_tm) uwerr(phi4_tm[i],wpm) end
			for i in 1:length(fpi_tm) uwerr(fpi_tm[i],wpm) end
			for i in 1:length(fk_tm) uwerr(fk_tm[i],wpm) end

			#=========================== FVE correction ==========================================================================================#

				m = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
				n = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

				if switch_FVE == 1
					for j in 1:length(mpi_tm)
						ji_pi = value((mpi_tm[j] / (4 * pi * fpi_tm[j])) ^2)
						mpiL = value(mpi_tm[j] * L[index])
						ji_k = value((mk_tm[j] / (4 * pi * fk_tm[j])) ^2)
						mkL = value(mk_tm[j] * L[index])
						lam_pi = [mpiL * sqrt(n[i]) for i in 1:length(n)]
						lam_k = [mkL * sqrt(n[i]) for i in 1:length(n)]
						g1_vec_pi = 4 .* m ./ lam_pi .* besselk.(1,lam_pi)
						g1_pi = sum(g1_vec_pi)
						g1_vec_k = 4 .* m ./ lam_k .* besselk.(1,lam_k)
						g1_k = sum(g1_vec_k)
				
						FVE_mpi =  2 * (0.25 * ji_pi * g1_pi)
						FVE_fpi =  2 * (-ji_pi * g1_pi - 1 / 2 * ji_k * g1_k)
						mpi_tm[j] = mpi_tm[j] / (1 + FVE_mpi)
						fpi_tm[j] = fpi_tm[j] / (1 + FVE_fpi)
				
						FVE_fk =  2 * (-3 / 8 * ji_pi * g1_pi - 3 / 4 * ji_k * g1_k)
						fk_tm[j] = fk_tm[j] / (1 + FVE_fk)
					end	
				end
				for i in 1:length(mpi_tm) uwerr(mpi_tm[i],wpm) end
				for i in 1:length(fpi_tm) uwerr(fpi_tm[i],wpm) end
				for i in 1:length(fk_tm) uwerr(fk_tm[i],wpm) end

                mpi_tm_sh = [mpi_tm[i] + dm * ( der[2] ) for i in 1:length(mpi_tm)]
			syst_mpi_sh = [syst_mpi[i] + dm * ( der[12] ) for i in 1:length(syst_mpi)]
			mk_tm_sh = [mk_tm[i] + dm * ( der[4] ) for i in 1:length(mk_tm)]
			syst_mk_sh = [syst_mk[i] + dm * ( der[14] ) for i in 1:length(syst_mk)]
			m12_tm_sh = [m12_tm[i] + dm * ( der[6] ) for i in 1:length(m12_tm)]
			syst_m12_sh = [syst_m12[i] + dm * ( der[16] ) for i in 1:length(syst_m12)]
			m12_tm_sh = [m12_tm_sh[i] for i in 1:length(m12_tm_sh)]#+ uwreal([0.0,value(syst_m12_sh[i])],string("systematic error m12_tm_",i," ",identity)) for i in 1:length(m12_tm_sh)]
			m13_tm_sh = [m13_tm[i] + dm * ( der[22] ) for i in 1:length(m13_tm)]
			syst_m13_sh = [syst_m13[i] + dm * ( der[26] ) for i in 1:length(syst_m13)]
			m13_tm_sh = [m13_tm_sh[i] for i in 1:length(m13_tm_sh)]#+ uwreal([0.0,value(syst_m13_sh[i])],string("systematic error m13_tm_",i," ",identity)) for i in 1:length(m13_tm_sh)]
			m34_tm_sh = [m34_tm[i] + dm * ( der[24] ) for i in 1:length(m34_tm)]
			syst_m34_sh = [syst_m34[i] + dm * ( der[28] ) for i in 1:length(syst_m34)]
			m34_tm_sh = [m34_tm_sh[i] for i in 1:length(m34_tm_sh)]#+ uwreal([0.0,value(syst_m34_sh[i])],string("systematic error m34_tm_",i," ",identity)) for i in 1:length(m34_tm_sh)]
			fpi_tm_sh = [fpi_tm[i] + dm * ( der[8] ) for i in 1:length(fpi_tm)]
			Rpi_tm_sh = fpi_tm_sh .* mpi_tm_sh .^ 1.5 ./ (sqrt(2) .* (mu1_pi .+ mu2_pi))	
			syst_Rpi_sh = [syst_Rpi[i] + dm * ( der[18] ) for i in 1:length(syst_Rpi)]
			Rpi_tm_sh = [Rpi_tm_sh[i] for i in 1:length(Rpi_tm_sh)]#+ uwreal([0.0,value(syst_Rpi_sh[i])],string("systematic error Rpi_tm_",i," ",identity)) for i in 1:length(Rpi_tm_sh)]
			fk_tm_sh = [fk_tm[i] + dm * ( der[10]) for i in 1:length(fk_tm)]
			Rk_tm_sh = fk_tm_sh .* mk_tm_sh .^ 1.5 ./ (sqrt(2) .* (mu1_k .+ mu2_k))	
			syst_Rk_sh = [syst_Rk[i] + dm * ( der[20] ) for i in 1:length(syst_Rk)]
			Rk_tm_sh = [Rk_tm_sh[i] for i in 1:length(Rk_tm_sh)]#+ uwreal([0.0,value(syst_Rk_sh[i])],string("systematic error Rk_tm_",i," ",identity)) for i in 1:length(Rk_tm_sh)]
			mpi_tm_sh = [mpi_tm_sh[i] for i in 1:length(mpi_tm_sh)]#+ uwreal([0.0,value(syst_mpi_sh[i])],string("systematic error mpi_tm_",i," ",identity)) for i in 1:length(mpi_tm_sh)]
			mk_tm_sh = [mk_tm_sh[i] for i in 1:length(mk_tm_sh)]#+ uwreal([0.0,value(syst_mk_sh[i])],string("systematic error mk_tm_",i," ",identity)) for i in 1:length(mk_tm_sh)]
			fpi_tm_sh = Rpi_tm_sh .* sqrt(2) .* (mu1_pi .+ mu2_pi) ./ mpi_tm_sh .^ 1.5
			fk_tm_sh = Rk_tm_sh .* sqrt(2) .* (mu1_k .+ mu2_k) ./ mk_tm_sh .^ 1.5

			phi2_tm_sh = [8*t0_sh[index] *mpi_tm_sh[i] ^2 for i in 1:length(mpi_tm_sh)]
			phi4_tm_sh = Array{uwreal,1}(undef, length(mk_tm_sh))

            count = 0
			if identity == "D200" || identity == "J303"
				for i in 1:length(mpi_tm_sh)
					for j in i+count:i+count+1
						phi4_tm_sh[j] = 8*t0_sh[index] *(0.5 *mpi_tm_sh[i] ^2 + mk_tm_sh[j] ^2)
					end
					count = count + 1
				end
			else
				for i in 1:length(mpi_tm_sh)
					for j in i+count:i+count+2
						phi4_tm_sh[j] = 8*t0_sh[index] *(0.5 *mpi_tm_sh[i] ^2 + mk_tm_sh[j] ^2)
					end
					count = count + 2
				end
			end
			
		
			for i in 1:length(mpi_tm_sh) uwerr(mpi_tm_sh[i],wpm) end
			for i in 1:length(mk_tm_sh) uwerr(mk_tm_sh[i],wpm) end
			for i in 1:length(m12_tm_sh) uwerr(m12_tm_sh[i],wpm) end
			for i in 1:length(m13_tm_sh) uwerr(m13_tm_sh[i],wpm) end
			for i in 1:length(m34_tm_sh) uwerr(m34_tm_sh[i],wpm) end
			for i in 1:length(phi2_tm_sh) uwerr(phi2_tm_sh[i],wpm) end
			for i in 1:length(phi4_tm_sh) uwerr(phi4_tm_sh[i],wpm) end
			for i in 1:length(fpi_tm_sh) uwerr(fpi_tm_sh[i],wpm) end
			for i in 1:length(fk_tm_sh) uwerr(fk_tm_sh[i],wpm) end




            y = [m12_tm_sh; phi2_tm_sh; phi4_tm_sh]
			target_m12 = [uwreal([0.0,0.0], "m12_target") for i in 1:length(m12_tm_sh)]
			target_phi2 = [phi2_sh[INDEX] for i in 1:length(phi2_tm_sh)]
			target_phi4 = [phi4_ph for i in 1:length(phi4_tm_sh)]
			target = [target_m12; target_phi2; target_phi4]
			y = y .- target
			for i in 1:length(y) uwerr(y[i],wpm) end
			dy = err.(y)
			W = 1 ./dy.^2
			
			println("performing matching fit...")
			count = 1
			count2 = 1

            p00 = rand(11)
							p00[9] = (minimum(x1[:,1])+maximum(x1[:,1]))/2
							if ens_tm[INDEX] == "H105r005_tm"
								p00[9] = 0.13730834834687144
							end 
							p00[10] = mul_guess[INDEX]
							p00[11] = mus_guess[INDEX]
							lb = [-Inf for i in 1:length(p00)]
							ub = [+Inf for i in 1:length(p00)]
							lb[9] = minimum(x1[:,1]) - 0.001
							ub[9] = maximum(x1[:,1]) + 0.001
							lb[10] = minimum(x1[:,2]) - 0.0001
							ub[10] = maximum(x1[:,2]) + 0.0001
							lb[11] = minimum(x3[:,3]) - 0.001 
							ub[11] = maximum(x3[:,3]) + 0.001 

							global chisq = fit_defs(model2,x,W)
							global fit = curve_fit(model2,x,value.(y),W,p00,lower=lb,upper=ub)
							
                            (up,chi_exp) = fit_error(chisq,coef(fit),y,wpm)

                            println("fit_error for the match accomplished")
			for i in 1:length(up) uwerr(up[i],wpm) end #uwerr.(up) 
			#details.([up[8],up[9],up[12]])
			#println("chisq/chiexp for matching:	",sum(fit.resid.^2)," / ",chi_exp," = ", sum(fit.resid.^2) / chi_exp ," (d.o.f.: ",dof(fit),")")
			Xexp = sum(fit.resid.^2) / chi_exp

			mul = up[10]
			muls = ( up[10] + up[11] ) / 2
			push!(kappa_c, up[9])
			uwerr.(kappa_c)
			uwerr(mul)
			uwerr(muls)

			#======================================= extrapolate fpi, fk to matching point ==========================================================#
				
			aux_mul = x1[:,2] 
			if ens[INDEX] == "D200" || ens[INDEX] == "J303"
				aux_mus = [x3[1:2,3]; x3[5:6,3]; x3[9:10,3]]
			else
				aux_mus = [x3[1:3,3]; x3[10:12,3]; x3[19:21,3]]
			end
			aux_kl = x1[:,1] 


            println("extrapolating for fpi...", Dates.format(now(), "yyyy-mm-dd at HH:MM"))	    
			count = 1
			warning = 1
			while count != 0
				try 
					p00 = rand(4)

					if ens[INDEX] == "N200"
						p00 = [-604.3692435328087, 1.675362269196541, 166.3713254709107, -11.448862401828348] .+ rand(d,4)
						#[-589.6948052900932, 1.4925579748368099, 162.33115926911583, -11.1707718050524]
					end 

					if ens[INDEX] == "D200"
						p00 = [-16353.517022555077, 3.9567291603839196, 4492.010313762788, -308.46739315325317] .+ rand(d,4)
					end 

					if ens[INDEX] == "N203"
						p00 = [-7326.663413882019, 1.2163284271786738, 2012.2829557409045, -138.16856418510395] .+ rand(d,4)
					end 

					if ens[INDEX] == "J303"
						p00 = [-36491.886422905816, 1.3770688409214304, 10015.502866417375, -687.2089564965605] .+ rand(d,4)
					end 

					global chisq = fit_defs(model_fpi, [aux_kl aux_mul], 1 ./ err.(fpi_tm_sh) .^2)
					global fit = curve_fit(model_fpi, [aux_kl aux_mul], value.(fpi_tm_sh), 1 ./ err.(fpi_tm_sh) .^2, p00)
						
					count = 0
				catch e
					if warning == 1
						@warn string("bad p0 guess in ", ens[INDEX], " fpi extrapolation")
					end
					warning = 0
				end
			end
			println("extrapolation for fpi done, now propagating errors...", Dates.format(now(), "yyyy-mm-dd at HH:MM"))

			if ens[INDEX] == "J303"
				wpm["J303"] = [-1.0, -1.0, 3.0, -1.0]
				(fpi_parm,chi_exp) = fit_error(chisq,coef(fit),fpi_tm_sh,wpm)
			else
				(fpi_parm,chi_exp) = fit_error(chisq,coef(fit),fpi_tm_sh,wpm)
			end
			println("fit_error done")
			Xexp_pi = sum(fit.resid.^2) / chi_exp    
			p = fpi_parm
			println("fpi_param = ", value.(fpi_parm))
			FPI = p[1] + p[2] * up[10] + p[3] / up[9] + p[4] / up[9]^2 
			uwerr(FPI,wpm)
			
			#=============== fk extrapolation ======================#
			
			println("extrapolating fk...", Dates.format(now(), "yyyy-mm-dd at HH:MM"))	    
			count = 1
			warning = 1
			while count != 0
				try 
					p00 = rand(5)

					if ens[INDEX] == "N200"
						p00 = [-5095.405652121227, 0.7051131763545394, 0.4643100248473917, 1399.4613925868366, -96.09024695413206] .+ rand(d,5)
					end

					if ens[INDEX] == "D200"
						p00 = [-5868.559627493074, 1.333155740521473, 0.4477190112564275, 1611.6841682131824, -110.65349447807347] .+ rand(d,5)
					end

					if ens[INDEX] == "N203"
						p00 = [-3954.3002256675873, 0.6238365048162324, 0.5087213437571475, 1086.0621292730473, -74.57183107512363] .+ rand(d,5)
					end

					if ens[INDEX] == "J303"
						p00 = [-4909.9742799954665, 0.41398865635161347, 0.2851122040210316, 1347.7873836412173, -92.49118556713438] .+ rand(d,5)
					end 

					global chisq = fit_defs(model_fk, x3, 1 ./ err.(fk_tm_sh) .^2)
					global fit = curve_fit(model_fk, x3, value.(fk_tm_sh), 1 ./ err.(fk_tm_sh) .^2, p00)
						
					count = 0
				catch e
					if warning == 1
						@warn string("bad p0 guess in ", ens[INDEX], " fk extrapolation")
					end
					warning = 0
				end
			end
			println("extrapolation fit for fk done, now propagating errors...", Dates.format(now(), "yyyy-mm-dd at HH:MM"))
					
			(fk_parm,chi_exp) = fit_error(chisq,coef(fit),fk_tm_sh,wpm)
			println("fit_error done")
			Xexp_k = sum(fit.resid.^2) / chi_exp
			q = fk_parm
			println("fk_param = ", value.(fk_parm))
			FK = q[1] + q[2] * up[10] + q[3] * up[11] + q[4] / up[9] + q[5] / up[9]^2 
			uwerr(FK,wpm)	



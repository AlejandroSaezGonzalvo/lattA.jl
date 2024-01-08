using juobs

function get_corr_wil(path::String, ens::EnsInfo, g1::String, g2::String; rw=false, info=false, legacy=false, fs=false)
    path_rw = joinpath(path, ens.id, "rwf")
    path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    path = joinpath(path, ens.id, "wil")
    path = filter(x->occursin(".mesons.dat", x), readdir(path, join=true))

    rwf = read_ms1.(path_rw, v=ens.vrw)
    dat = read_mesons([path[i] for i in 1:length(path)], g1, g2, legacy=legacy)
    truncate_data!(dat,ens.cnfg)

    rw ? corr = [corr_obs(dat[i], L=ens.L, rw=rwf, info=info) for i in 1:length(dat)] : corr = [corr_obs(dat[i], L=L[index], info=info, flag_strange=fs) for i in 1:length(dat)]

    if info == false
        return corr
    else
        pp = getindex.(corr,1)
        ppw = getindex.(corr,2)
        w = getindex.(corr,3)
        return pp, ppw, w[1]
    end
end

function get_corr_tm(path::String, ens::EnsInfo, g1::String, g2::String; rw=false, info=false, legacy=false, fs=false)
    path_rw = joinpath(path, ens.id, "rwf")
    path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    path = joinpath(path, ens.id, "tm")
    path = filter(x->occursin(".mesons.dat", x), readdir(path, join=true))

    rwf = read_ms1.(path_rw, v=ens.vrw)
    dat = read_mesons([path[i] for i in 1:length(path)], g1, g2, legacy=legacy)
    truncate_data!(dat,ens.cnfg)

    rw ? corr = [corr_obs(dat[i], L=ens.L, rw=rwf, info=info) for i in 1:length(dat)] : corr = [corr_obs(dat[i], L=L[index], info=info, flag_strange=fs) for i in 1:length(dat)]

    if info == false
        return corr
    else
        pp = getindex.(corr,1)
        ppw = getindex.(corr,2)
        w = getindex.(corr,3)
        return pp, ppw, w[1]
    end
end

function get_dSdm(path::String, ens::EnsInfo)
    path = joinpath(path, ens.id, "dSdm")
    path = filter(x->occursin(".dat", x), readdir(path, join=true))
    dSdm = read_md.(path)
    return dSdm
end

function read_ens_wil(path::String, ens::EnsInfo; legacy=false, fs=false)
    pp, ppw, w = get_corr_wil(path, ens, "G5", "G5", rw=true, info=true, legacy=legacy);
    pp_sym = [corr_sym(pp[i], pp[i+1], +1) for i in 1:2:length(pp)-1];
    ap, apw, w = get_corr_wil(path, ens, "G5", "G0G5", rw=true, info=true, legacy=legacy);
    ap_sym = [corr_sym(ap[i], ap[i+1], -1) for i in 1:2:length(ap)-1];

    pp_d1 = get_corr_wil(path, ens, "G5_d1", "G5_d1", rw=true, legacy=legacy);
    pp_d2 = get_corr_wil(path, ens, "G5_d2", "G5_d2", rw=true, legacy=legacy);
    ap_d1 = get_corr_wil(path, ens, "G5_d1", "G0G5_d1", rw=true, legacy=legacy);
    ap_d2 = get_corr_wil(path, ens, "G5_d2", "G0G5_d2", rw=true, legacy=legacy);
    dSdm = get_dSdm(path, ens)

    pp_val = [[pp_d1[i], pp_d2[i]] for i in 1:length(pp_d1)];
    ap_val = [[ap_d1[i], ap_d2[i]] for i in 1:length(ap_d1)];
    corr = [[pp[i] for i in 1:length(pp)]; [ap[i] for i in 1:length(ap)]];
    corr_val = [[pp_val[i] for i in 1:length(pp)]; [ap_val[i] for i in 1:length(ap)]];
    corrw = [[ppw[i] for i in 1:length(pp)]; [apw[i] for i in 1:length(ap)]];

    return pp_sym, ap_sym, corr, corr_val, corrw, dSdm
end

function read_ens_tm_sym(path::String, ens::EnsInfo; legacy=false)
    pp, ppw, w = get_corr_tm(path, ens, "G5", "G5", rw=true, info=true, legacy=legacy);
    pp_sym = [corr_sym(pp[i], pp[i+9], +1) for i in 1:9];
    ap, apw, w = get_corr_tm(path, ens, "G5", "G0G5", rw=true, info=true, legacy=legacy);
    ap_sym = [corr_sym(ap[i], ap[i+9], -1) for i in 1:9];

    dSdm = get_dSdm(path, ens)

    corrw = [[ppw[i] for i in 1:length(pp)]; [apw[i] for i in 1:length(ap)]];

    return pp_sym, ap_sym, corrw, dSdm
end

function read_ens_tm(path::String, ens::EnsInfo; legacy=false)
    pp, ppw, w = get_corr_tm(path, ens, "G5", "G5", rw=true, info=true, legacy=legacy);
    pp_sym = [corr_sym(pp[i], pp[i+24], +1) for i in 1:24];
    ap, apw, w = get_corr_tm(path, ens, "G5", "G0G5", rw=true, info=true, legacy=legacy);
    ap_sym = [corr_sym(ap[i], ap[i+24], -1) for i in 1:24];

    dSdm = get_dSdm(path, ens)

    corrw = [[ppw[i] for i in 1:length(pp)]; [apw[i] for i in 1:length(ap)]];

    return pp_sym, ap_sym, corrw, dSdm
end

function read_ens_csv(path::String, ens::EnsInfo)
    ix = ensemble_inv[ens.id]
    path_ll = Path_ll[ix]
    path_ls = Path_ls[ix]
    path_ss = Path_ss[ix]
    path_ll_sym = Path_ll_sym[ix]
    path_ls_sym = Path_ls_sym[ix]
    path_ss_sym = Path_ss_sym[ix]
    path2_ll = Path2_ll[ix]
    path2_ls = Path2_ls[ix]
    path2_ss = Path2_ss[ix]
    path2_ll_sym = Path2_ll_sym[ix]
    path2_ls_sym = Path2_ls_sym[ix]
    path2_ss_sym = Path2_ss_sym[ix]

    pp_ll = [csv2Corr(path_ll[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path_ll)]
	pp_ls = [csv2Corr(path_ls[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path_ls)]
    pp_ss = [csv2Corr(path_ss[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path_ss)]
    ap_ll = [csv2Corr(path2_ll[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path2_ll)] 
	ap_ls = [csv2Corr(path2_ls[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path2_ls)] 
    ap_ss = [csv2Corr(path2_ss[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path2_ss)]

    pp_ll_2 = [csv2Corr(path_ll_sym[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path_ll_sym)] 
	pp_ls_2 = [csv2Corr(path_ls_sym[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path_ls_sym)]
    pp_ss_2 = [csv2Corr(path_ss_sym[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path_ss_sym)]
    ap_ll_2 = [csv2Corr(path2_ll_sym[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path2_ll_sym)] 
	ap_ls_2 = [csv2Corr(path2_ls_sym[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path2_ls_sym)]
	ap_ss_2 = [csv2Corr(path2_ss_sym[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path2_ss_sym)]

	pp_ll_sym = [corr_sym(pp_ll[i],pp_ll_2[i],+1) for i in 1:1:length(pp_ll)]
	pp_ls_sym = [corr_sym(pp_ls[i],pp_ls_2[i],+1) for i in 1:1:length(pp_ls)]
	ap_ll_sym = [corr_sym(ap_ll[i],ap_ll_2[i],-1) for i in 1:1:length(ap_ll)]
	ap_ls_sym = [corr_sym(ap_ls[i],ap_ls_2[i],-1) for i in 1:1:length(ap_ls)]
	pp_ss_sym = [corr_sym(pp_ss[i],pp_ss_2[i],+1) for i in 1:1:length(pp_ss)]
	ap_ss_sym = [corr_sym(ap_ss[i],ap_ss_2[i],-1) for i in 1:1:length(ap_ss)]

    i=j=1
    pp_sym = Array{juobs.Corr,1}()
    ap_sym = Array{juobs.Corr,1}()
    while i < length(pp_ll_sym)
        pp_sym = vcat(pp_sym, [pp_ll_sym[i:i+2]; pp_ls_sym[j:j+8]; pp_ss_sym[i:i+2]])
        ap_sym = vcat(ap_sym, [ap_ll_sym[i:i+2]; ap_ls_sym[j:j+8]; ap_ss_sym[i:i+2]])
        i+=3
        j+=9
    end

    return pp_sym, ap_sym
end

function get_YW(path::String, ens::EnsInfo, plat::Vector{Int64}; rw=false, npol::Int64=2, ws::ADerrors.wspace=ADerrors.wsg)

    println("WARNING!: You must use the same plat here than the one used to compute t0")

    path_ms = joinpath(path, ens.id, "gf")
    path_ms = filter(x->occursin(".dat", x), readdir(path_ms, join=true))
    Y = read_ms.(path_ms, dtr=ens.dtr) 

    nr = length(Y)
    Ysl = getfield.(Y, :obs)
    t = getfield.(Y, :t)
    t = t[1]
    id = getfield.(Y, :id)
    replica = size.(Ysl, 1)

    L = ens.L
    id = ens.id
    #T = length(Y[:,1]) - y0
    y0 = 1 ## assumes this is the case, hardcoded, some ensembles will not fulfil !
    println("WARNING!: make sure t_src is 1 in this ensemble")

    #Truncation
    if id in keys(ADerrors.wsg.str2id)
        n_ws = findfirst(x-> x == ws.str2id[id], ws.map_nob)
        if !isnothing(n_ws)
            ivrep_ws = ws.fluc[n_ws].ivrep

            if length(replica) != length(ivrep_ws)
                error("Different number of replicas")
            end

            for k = 1:length(replica)
                if replica[k] > ivrep_ws[k]
                    println("Automatic truncation in Ysl ", ivrep_ws[k], " / ", replica[k], ". R = ", k)
                    Ysl[k] = Ysl[k][1:ivrep_ws[k], :, :]
                elseif replica[k] < ivrep_ws[k]
                    error("Automatic truncation failed. R = ", replica[k], "\nTry using truncate_data!")
                end
            end
            replica = size.(Ysl, 1)
        end
    end
    
    tmp = Ysl[1]
    [tmp = cat(tmp, Ysl[k], dims=1) for k = 2:nr]
    nt0 = juobs.t0_guess(t, tmp, plat, L)
    xmax = size(tmp, 2)
    T = xmax - 1 - y0

    dt0 = iseven(npol) ? Int64(npol / 2) : Int64((npol+1) / 2)
    Y_aux = Matrix{uwreal}(undef, xmax, 2*dt0+1)

    if rw
        path_rw = joinpath(path, ens.id, "rwf")
        path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
        rwf = read_ms1.(path_rw)

        Ysl_r, W = juobs.apply_rw(Ysl, rwf)
        tmp_r = Ysl_r[1]
        tmp_W = W[1]
        [tmp_r = cat(tmp_r, Ysl_r[k], dims=1) for k = 2:nr]
        [tmp_W = cat(tmp_W, W[k], dims=1) for k = 2:nr]
        W_obs = uwreal(tmp_W, id, replica)
        WY_aux = Matrix{uwreal}(undef, xmax, 2*dt0+1)
    end
    for i = 1:xmax
        k = 1
        for j = nt0-dt0:nt0+dt0
            if !rw
                Y_aux[i, k] = uwreal(tmp[:, i, j], id, replica)
            else
                WY_aux[i, k] = uwreal(tmp_r[:, i, j], id, replica)
                Y_aux[i, k] = WY_aux[i, k] / W_obs
            end
            k = k + 1
        end
    end

    return WY_aux, W_obs
end

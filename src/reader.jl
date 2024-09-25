using juobs

include("/home/asaez/cls_ens/codes/lattA.jl/src/path_csv.jl");


function read_mesons_multichunks(path::String, ens::String, g1::String="G5", g2::String="G5"; legacy::Bool=false)
    idx_ts001 = findall(x->occursin("ts001",x), db[ens])
    idx_tsT = findall(x->!occursin("ts001",x), db[ens])

    store_cdata_aux = []
    for i in db[ens]
        aux = filter(x->occursin(i, x), readdir(path, join=true))
        data_chunk = read_mesons(aux, g1, g2, legacy=legacy);
        push!(store_cdata_aux, data_chunk)
    end

    #=
    if ens == "E300"
        for i in 1:length(idx_ts001)-1
            concat_data!(store_cdata_aux[idx_ts001[1]][1:2:end], store_cdata_aux[idx_ts001[i+1]])    
        end
        concat_data!(store_cdata_aux[idx_tsT[1]], store_cdata_aux[idx_ts001[1]][2:2:end]) 
        for i in 1:length(idx_tsT)-1
            concat_data!(store_cdata_aux[idx_tsT[1]], store_cdata_aux[idx_tsT[i+1]])    
        end
        store_cdata_aux[idx_ts001[1]] = store_cdata_aux[idx_ts001[1]][1:2:end]
    else
        for i in 1:length(idx_ts001)-1
            concat_data!(store_cdata_aux[idx_ts001[1]], store_cdata_aux[idx_ts001[i+1]])    
        end
        for i in 1:length(idx_tsT)-1
            concat_data!(store_cdata_aux[idx_tsT[1]], store_cdata_aux[idx_tsT[i+1]])    
        end
    end
    =#
    for i in 1:length(idx_ts001)-1
        concat_data!(store_cdata_aux[idx_ts001[1]], store_cdata_aux[idx_ts001[i+1]])    
    end
    for i in 1:length(idx_tsT)-1
        concat_data!(store_cdata_aux[idx_tsT[1]], store_cdata_aux[idx_tsT[i+1]])    
    end
    dat_ts001 = store_cdata_aux[idx_ts001[1]]
    dat_ts190 = store_cdata_aux[idx_tsT[1]]
    
    dat = Vector{Vector{juobs.CData}}()
    for i in 1:length(dat_ts001)
        push!(dat, dat_ts001[i])
        push!(dat, dat_ts190[i])
    end

    return dat
end

function read_mesons_correction_multichunks(path::String, ens::String, g1::String="G5", g2::String="G5"; legacy::Bool=false)
    idx_ts001 = findall(x->occursin("ts001",x), db_c[ens])
    idx_tsT = findall(x->!occursin("ts001",x), db_c[ens])

    store_cdata_aux = []
    for i in db_c[ens]
        aux = filter(x->occursin(i, x), readdir(path, join=true))
        data_chunk = read_mesons_correction(aux, g1, g2, legacy=legacy);
        push!(store_cdata_aux, data_chunk)
    end

    for i in 1:length(idx_ts001)-1
        concat_data!(store_cdata_aux[idx_ts001[1]], store_cdata_aux[idx_ts001[i+1]])    
    end
    for i in 1:length(idx_tsT)-1
        concat_data!(store_cdata_aux[idx_tsT[1]], store_cdata_aux[idx_tsT[i+1]])    
    end
    dat_ts001 = store_cdata_aux[idx_ts001[1]]
    if length(idx_tsT) >= 1
        dat_ts190 = store_cdata_aux[idx_tsT[1]]
        dat = Vector{Vector{juobs.CData}}()
        for i in 1:length(dat_ts001)
            push!(dat, dat_ts001[i])
            push!(dat, dat_ts190[i])
        end
    else
        dat = dat_ts001
    end

    return dat
end

function get_corr_TSM_multichunks(path::String, ens::EnsInfo; info=false)
    path = joinpath(path, ens.id)
    path_sl = joinpath.(path, "sloppy")
    path_c = joinpath.(path, "correc")
    path_rw = joinpath(path, "rwf_def")
    path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    if length(path_rw) == 0
        path_rw = joinpath(path, "rwf")
        global path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    end
    
    pp_dat = read_mesons_multichunks(path_sl, ens.id, "G5", "G5")
    ap_dat = read_mesons_multichunks(path_sl, ens.id, "G5", "G0G5")
    pp_dat_c = read_mesons_correction_multichunks(path_c, ens.id, "G5", "G5")
    ap_dat_c = read_mesons_correction_multichunks(path_c, ens.id, "G5", "G0G5")

    rwf = [read_ms1(path_rw[i], v=ens.vrw[i]) for i in 1:length(ens.vrw)]
    cnfg_rw = size.(rwf,2)
    cnfg_trunc_ts001 = [findall(pp_dat[1][i].vcfg .< cnfg_rw[i])[end] for i in 1:length(cnfg_rw)] ## rwf missing some configs at the end
    cnfg_trunc_ts001_c = [findall(pp_dat_c[1][i].vcfg .< cnfg_rw[i])[end] for i in 1:length(cnfg_rw)]
    cnfg_trunc_ts190 = [findall(pp_dat[2][i].vcfg .< cnfg_rw[i])[end] for i in 1:length(cnfg_rw)] ## rwf missing some configs at the end
    cnfg_trunc_ts190_c = [findall(pp_dat_c[2][i].vcfg .< cnfg_rw[i])[end] for i in 1:length(cnfg_rw)]
    truncate_data!(pp_dat[1:2:end], cnfg_trunc_ts001)
    truncate_data!(ap_dat[1:2:end], cnfg_trunc_ts001)
    truncate_data!(pp_dat_c[1:2:end], cnfg_trunc_ts001_c)
    truncate_data!(ap_dat_c[1:2:end], cnfg_trunc_ts001_c)
    truncate_data!(pp_dat[2:2:end], cnfg_trunc_ts190) 
    truncate_data!(ap_dat[2:2:end], cnfg_trunc_ts190)
    truncate_data!(pp_dat_c[2:2:end], cnfg_trunc_ts190_c)
    truncate_data!(ap_dat_c[2:2:end], cnfg_trunc_ts190_c)

    if sym_bool[ens.id] == true
        pp = corr_obs_TSM.(pp_dat[1:length(ap_dat_c)], pp_dat_c[1:length(ap_dat_c)], rw=rwf, L=ens.L, info=info, replica_sl=ens.cnfg, nms=sum(ens.cnfg))
        ap = corr_obs_TSM.(ap_dat[1:length(ap_dat_c)], ap_dat_c[1:length(ap_dat_c)], rw=rwf, L=ens.L, info=info, replica_sl=ens.cnfg, nms=sum(ens.cnfg))
    
        if ens.id == "E250"
            pp_sym = [corr_sym_E250(pp[i], pp[i+1], +1) for i in 1:2:length(ap)]
            ap_sym = [corr_sym_E250(ap[i], ap[i+1], -1) for i in 1:2:length(ap)]
        elseif ens.id == "D450"
            pp_sym = [corr_sym_D450(pp_ts001[i], pp_tsT[i], +1) for i in 1:length(pp_ts001)]
            ap_sym = [corr_sym_D450(ap_ts001[i], ap_tsT[i], -1) for i in 1:length(pp_ts001)]
        else
            pp_sym = [corr_sym(pp[i], pp[i+1], +1) for i in 1:2:length(ap)]
            ap_sym = [corr_sym(ap[i], ap[i+1], -1) for i in 1:2:length(ap)]
        end
    else
        pp_ts001 = corr_obs_TSM.(pp_dat[1:2:length(ap_dat)], pp_dat_c[1:length(ap_dat_c)], rw=rwf, L=ens.L, info=info, replica_sl=ens.cnfg, nms=sum(ens.cnfg))
        ap_ts001 = corr_obs_TSM.(ap_dat[1:2:length(ap_dat)], ap_dat_c[1:length(ap_dat_c)], rw=rwf, L=ens.L, info=info, replica_sl=ens.cnfg, nms=sum(ens.cnfg))
        pp_tsT = corr_obs.(pp_dat[2:2:length(ap_dat)], rw=rwf, L=ens.L, info=info, replica=ens.cnfg, nms=sum(ens.cnfg))
        ap_tsT = corr_obs.(ap_dat[2:2:length(ap_dat)], rw=rwf, L=ens.L, info=info, replica=ens.cnfg, nms=sum(ens.cnfg))
    
        if ens.id == "E250"
            pp_sym = [corr_sym_E250(pp_ts001[i], pp_tsT[i], +1) for i in 1:length(pp_ts001)]
            ap_sym = [corr_sym_E250(ap_ts001[i], ap_tsT[i], -1) for i in 1:length(pp_ts001)]
        elseif ens.id == "D450"
            pp_sym = [corr_sym_D450(pp_ts001[i], pp_tsT[i], +1) for i in 1:length(pp_ts001)]
            ap_sym = [corr_sym_D450(ap_ts001[i], ap_tsT[i], -1) for i in 1:length(pp_ts001)]
        else
            pp_sym = [corr_sym(pp_ts001[i], pp_tsT[i], +1) for i in 1:length(pp_ts001)]
            ap_sym = [corr_sym(ap_ts001[i], ap_tsT[i], -1) for i in 1:length(pp_ts001)]
        end
    end
    
    return pp_sym, ap_sym
end

function get_corr_TSM(path::String, ens::EnsInfo, g1::String, g2::String; rw=false, info=false, legacy=false, fs=false)
    path_rw = joinpath(path, ens.id, "rwf_def")
    path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    if length(path_rw) == 0
        path_rw = joinpath(path, ens.id, "rwf")
        global path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    end
    path_sl = joinpath(path, ens.id, "sloppy")
    path_sl = filter(x->occursin(".mesons.dat", x), readdir(path_sl, join=true))
    path_c = joinpath(path, ens.id, "correc")
    path_c = filter(x->occursin(".mesons.dat", x), readdir(path_c, join=true))

    rwf = read_ms1.(path_rw, v=ens.vrw)
    dat_sl = read_mesons([path_sl[i] for i in 1:length(path_sl)], g1, g2, legacy=legacy, id=ens.id)
    dat_c = read_mesons_correction([path_c[i] for i in 1:length(path_c)], g1, g2, legacy=legacy, id=ens.id)

    rw ? corr = [corr_obs_TSM(dat_sl[i], dat_c[i], L=ens.L, rw=rwf, info=info) for i in 1:length(dat_sl)] : corr = [corr_obs_TSM(dat_sl[i], dat_c[i], L=ens.L, info=info) for i in 1:length(dat_sl)]

    if info == false
        return corr
    else
        pp = getindex.(corr,1)
        ppw = getindex.(corr,2)
        w = getindex.(corr,3)
        return pp, ppw, w[1]
    end
end

function get_corr_wil(path::String, ens::EnsInfo, g1::String, g2::String; rw=false, info=false, legacy=false, fs=false)
    path_rw = joinpath(path, ens.id, "rwf_def")
    path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    if length(path_rw) == 0
        path_rw = joinpath(path, ens.id, "rwf")
        global path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    end
    path = joinpath(path, ens.id, "wil")
    path = filter(x->occursin(".mesons.dat", x), readdir(path, join=true))

    if ens.id == "D200"
        #rwf_1 = read_ms1.([path_rw[1]], v=ens.vrw)
        #rwf_2 = read_ms1.([path_rw[2]], v=ens.vrw)
        #rwf = [hcat(rwf_1[1],rwf_2[1])]
        rwf = read_ms1.(path_rw, v=ens.vrw)
        dat = read_mesons([path[1]], g1, g2, legacy=legacy, id=ens.id)
        dat_2 = read_mesons([path[2]], g1, g2, legacy=legacy, id=ens.id)
        concat_data!(dat,dat_2)
    else
        rwf = read_ms1.(path_rw, v=ens.vrw)
        dat = read_mesons([path[i] for i in 1:length(path)], g1, g2, legacy=legacy, id=ens.id)
        truncate_data!(dat,ens.cnfg)
    end

    rw ? corr = [corr_obs(dat[i], L=ens.L, rw=rwf, info=info, flag_strange=fs) for i in 1:length(dat)] : corr = [corr_obs(dat[i], L=ens.L, info=info) for i in 1:length(dat)]

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
    path_rw = joinpath(path, ens.id, "rwf_def")
    path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    if length(path_rw) == 0
        path_rw = joinpath(path, ens.id, "rwf")
        global path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    end
    path = joinpath(path, ens.id, "tm")
    path = filter(x->occursin(".mesons.dat", x), readdir(path, join=true))

    if ens.id == "J303"
        rwf = read_ms1.(path_rw, v=ens.vrw)
        dat = read_mesons([path[1]], g1, g2, legacy=legacy, id=ens.id)
        dat_2 = read_mesons([path[2]], g1, g2, legacy=legacy, id=ens.id)
        concat_data!(dat,dat_2)
    elseif ens.id == "D200"
        #rwf_1 = read_ms1.([path_rw[1]], v=ens.vrw)
        #rwf_2 = read_ms1.([path_rw[2]], v=ens.vrw)
        #rwf = [hcat(rwf_2[1],rwf_1[1])]
        rwf = read_ms1.(path_rw, v=ens.vrw)
        dat = read_mesons([path[1]], g1, g2, legacy=legacy, id=ens.id)
        dat_2 = read_mesons([path[2]], g1, g2, legacy=legacy, id=ens.id)
        dat_3 = read_mesons([path[3]], g1, g2, legacy=legacy, id=ens.id)
        concat_data!(dat,dat_3)
        concat_data!(dat,dat_2)
    elseif ens.id == "N300"
        rwf = read_ms1.(path_rw, v=ens.vrw)
        dat_1 = read_mesons([path[1]], g1, g2, legacy=legacy, id=ens.id)
        dat = read_mesons([path[2]], g1, g2, legacy=legacy, id=ens.id)
        truncate_data!(dat,[199])
        concat_data!(dat,dat_1)
    else
        rwf = read_ms1.(path_rw, v=ens.vrw)
        dat = read_mesons([path[i] for i in 1:length(path)], g1, g2, legacy=legacy, id=ens.id)
        truncate_data!(dat,ens.cnfg)
    end

    rw ? corr = [corr_obs(dat[i], L=ens.L, rw=rwf, info=info) for i in 1:length(dat)] : corr = [corr_obs(dat[i], L=ens.L, info=info) for i in 1:length(dat)]

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

function read_ens_TSM(path::String, ens::EnsInfo; legacy=false, fs=false)
    pp = get_corr_TSM(path, ens, "G5", "G5", rw=true, info=false, legacy=legacy);
    ap = get_corr_TSM(path, ens, "G5", "G0G5", rw=true, info=false, legacy=legacy);
    idx_wil = findall([pp[i].mu == [.0,.0] for i in 1:length(pp)]);
    idx_tm = findall([pp[i].mu[1] != .0 for i in 1:length(pp)]);

    pp_sym = [corr_sym(pp[i], pp[i+1], +1) for i in idx_wil[1]:2:idx_wil[end]-1];
    ap_sym = [corr_sym(ap[i], ap[i+1], -1) for i in idx_wil[1]:2:idx_wil[end]-1];
    pp_tm_sym = [corr_sym(pp[i], pp[i+div(length(idx_tm),2)], +1) for i in idx_tm[1]:idx_tm[1]+div(length(idx_tm),2)-1];
    ap_tm_sym = [corr_sym(ap[i], ap[i+div(length(idx_tm),2)], -1) for i in idx_tm[1]:idx_tm[1]+div(length(idx_tm),2)-1];

    pp_sym = [pp_sym[1:3]; pp_tm_sym]
    ap_sym = [ap_sym[1:3]; ap_tm_sym]

    return pp_sym, ap_sym
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
    if ens.id == "D200"
        dSdm = get_dSdm(path, ens)
        aux = [hcat(dSdm[3], dSdm[1])]
        aux = [hcat(aux[1], dSdm[2])]
        dSdm = aux
    else
        dSdm = get_dSdm(path, ens)
    end

    pp_val = [[pp_d1[i], pp_d2[i]] for i in 1:length(pp_d1)];
    ap_val = [[ap_d1[i], ap_d2[i]] for i in 1:length(ap_d1)];
    corr = [[pp[i] for i in 1:length(pp)]; [ap[i] for i in 1:length(ap)]];
    corr_val = [[pp_val[i] for i in 1:length(pp)]; [ap_val[i] for i in 1:length(ap)]];
    corrw = [[ppw[i] for i in 1:length(pp)]; [apw[i] for i in 1:length(ap)]];

    return pp_sym, ap_sym, corr, corr_val, corrw, dSdm, w
end

function read_ens_tm_sym(path::String, ens::EnsInfo; legacy=false)
    pp, ppw, w = get_corr_tm(path, ens, "G5", "G5", rw=true, info=true, legacy=legacy);
    pp_sym = [corr_sym(pp[i], pp[i+9], +1) for i in 1:9];
    ap, apw, w = get_corr_tm(path, ens, "G5", "G0G5", rw=true, info=true, legacy=legacy);
    ap_sym = [corr_sym(ap[i], ap[i+9], -1) for i in 1:9];

    dSdm = get_dSdm(path, ens)

    corrw = [[ppw[i] for i in 1:length(pp)]; [apw[i] for i in 1:length(ap)]];

    return pp_sym, ap_sym, corrw, dSdm, w
end

function read_ens_tm(path::String, ens::EnsInfo; legacy=false)
    pp, ppw, w = get_corr_tm(path, ens, "G5", "G5", rw=true, info=true, legacy=legacy);
    pp_sym = [corr_sym(pp[i], pp[i+24], +1) for i in 1:24];
    ap, apw, w = get_corr_tm(path, ens, "G5", "G0G5", rw=true, info=true, legacy=legacy);
    ap_sym = [corr_sym(ap[i], ap[i+24], -1) for i in 1:24];

    if ens.id == "D200"
        dSdm = get_dSdm(path, ens)
        aux = [hcat(dSdm[3], dSdm[1])]
        aux = [hcat(aux[1], dSdm[2])]
        dSdm = aux
    else
        dSdm = get_dSdm(path, ens)
    end

    corrw = [[ppw[i] for i in 1:length(pp)]; [apw[i] for i in 1:length(ap)]];

    return pp_sym, ap_sym, corrw, dSdm, w
end

function read_ens_csv(ens::EnsInfo)
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

    return pp_sym, ap_sym, [pp_ll;pp_ls;pp_ss;pp_ll_2;pp_ls_2;pp_ss_2], [ap_ll;ap_ls;ap_ss;ap_ll_2;ap_ls_2;ap_ss_2]
end

function get_YM(path::String, ens::EnsInfo; rw=false, ws::ADerrors.wspace=ADerrors.wsg, w0::Union{Float64, Nothing}=nothing)

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
    xmax = size(tmp, 2)
    T = xmax - 1 - y0

    Y_aux = Matrix{uwreal}(undef, xmax, length(t))

    if rw
        path_rw = joinpath(path, ens.id, "rwf")
        path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
        if ens.id == "D200"
            rwf_1 = read_ms1.([path_rw[1]], v=ens.vrw)
            rwf_2 = read_ms1.([path_rw[2]], v=ens.vrw)
            rwf = [hcat(rwf_2[1],rwf_1[1])]
        elseif ens.id in ["E250", "E300", "J500", "J501", "D450"]
            rwf = [read_ms1(path_rw[i], v=ens.vrw[i]) for i in 1:length(ens.vrw)]
            [Ysl[k] = Ysl[k][1:size(rwf[k],2), :, :] for k in 1:length(Ysl)]
        else
            rwf = read_ms1.(path_rw, v=ens.vrw)
        end

        Ysl_r, W = juobs.apply_rw(Ysl, rwf, id=id)
        tmp_r = Ysl_r[1]
        tmp_W = W[1]
        [tmp_r = cat(tmp_r, Ysl_r[k], dims=1) for k = 2:nr]
        [tmp_W = cat(tmp_W, W[k], dims=1) for k = 2:nr]
        W_obs = uwreal(tmp_W, id, replica, collect(1:length(tmp_W)), sum(replica))
        WY_aux = Matrix{uwreal}(undef, xmax, length(t))
    end
    for i = 1:xmax
        k = 1
        for j = 1:length(t)
            if !rw
                Y_aux[i, k] = uwreal(tmp[:, i, j], id, replica)
            else
                WY_aux[i, k] = uwreal(tmp_r[:, i, j], id, replica, collect(1:length(tmp_W)), sum(replica))
                Y_aux[i, k] = WY_aux[i, k] / W_obs
            end
            k = k + 1
        end
    end

    t2YM = similar(Y_aux)
    tdt2YM = similar(Y_aux)
    for i in 1:length(Y_aux[:,1])
        t2YM[i,:] = Y_aux[i,:] .* t .^ 2 ./ L ^ 3
    end
    for i in 1:length(t2YM[:,1])
        if isnothing(w0)
            tdt2YM[i,2:end-1] = [(t2YM[i,j+1] - t2YM[i,j-1]) / (t[j+1] - t[j-1]) * t[j] for j in 2:length(t2YM[i,:])-1]
            tdt2YM[i,1] = tdt2YM[i,end] = t2YM[i,1]
        else
            ixm = findmin(abs.(t .- (w0-0.5)))[2]
            ixM = findmin(abs.(t[1:end-1] .- (w0+0.5)))[2]
            tdt2YM[i,1:end] .= t2YM[i,2]
            tdt2YM[i,ixm:ixM] = [(t2YM[i,j+1] - t2YM[i,j-1]) / (t[j+1] - t[j-1]) * t[j] for j in ixm:ixM]
        end
    end

    return t2YM, tdt2YM, W_obs, t
end

function get_YM_dYM(path::String, ens::EnsInfo; rw=false, ws::ADerrors.wspace=ADerrors.wsg, w0::Union{Float64, Nothing}=nothing, tau::Union{Float64, Nothing}=nothing)

    path_ms = joinpath(path, ens.id, "gf")
    path_ms = filter(x->occursin(".dat", x), readdir(path_ms, join=true))
    Y = read_ms.(path_ms, dtr=ens.dtr) 

    truncate_data!(Y, ens.cnfg)

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
    xmax = size(tmp, 2)
    T = xmax - 1 - y0

    Y_aux = Matrix{uwreal}(undef, xmax, length(t))

    if rw
        path_rw = joinpath(path, ens.id, "rwf")
        path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
        if ens.id in ["H105", "J500", "J501"]
            rwf = [read_ms1(path_rw[i], v=ens.vrw[i]) for i in 1:length(ens.vrw)]
            [Ysl[k] = Ysl[k][1:size(rwf[k],2), :, :] for k in 1:length(Ysl)]
        else
            rwf = read_ms1.(path_rw, v=ens.vrw)
        end

        Ysl_r, W = juobs.apply_rw(Ysl, rwf, id=id)
        tmp_r = Ysl_r[1]
        tmp_W = W[1]
        [tmp_r = cat(tmp_r, Ysl_r[k], dims=1) for k = 2:nr]
        [tmp_W = cat(tmp_W, W[k], dims=1) for k = 2:nr]
        W_obs = uwreal(tmp_W, id, replica, collect(1:length(tmp_W)), sum(replica))
        WY_aux = Matrix{uwreal}(undef, xmax, length(t))
    end
    for i = 1:xmax
        k = 1
        for j = 1:length(t)
            if !rw
                Y_aux[i, k] = uwreal(tmp[:, i, j], id, replica)
            else
                WY_aux[i, k] = uwreal(tmp_r[:, i, j], id, replica, collect(1:length(tmp_W)), sum(replica))
                Y_aux[i, k] = WY_aux[i, k] / W_obs
            end
            k = k + 1
        end
    end

    if tau == nothing
        t2YM = similar(Y_aux)
        tdt2YM = Matrix{uwreal}(undef,size(t2YM)[1],size(t2YM)[2]-2)
    else
        tau_ind = Int64(div(tau, t[2]-t[1]))
        if tau >= 0.0
            global t = t[1:end-tau_ind]
        elseif tau < 0.0
            global t = t[1-tau_ind:end]
        end
        t2YM = Matrix{uwreal}(undef,size(Y_aux)[1],length(t))
        tdt2YM = Matrix{uwreal}(undef,size(t2YM)[1],size(t2YM)[2]-2)
    end
    for i in 1:length(Y_aux[:,1])
        if tau == nothing
            t2YM[i,:] = Y_aux[i,:] .* t .^ 2 ./ L ^ 3
        else
            if tau >= 0.0
                t2YM[i,1:end] = Y_aux[i,1+tau_ind:end] .* t .^ 2 ./ L ^ 3
            elseif tau < 0.0
                t2YM[i,1:end] = Y_aux[i,1:end+tau_ind] .* t .^ 2 ./ L ^ 3
            end
        end
    end
    for i in 1:length(t2YM[:,1])
        if isnothing(w0)
            tdt2YM[i,1:end] = [(t2YM[i,j+1] - t2YM[i,j-1]) / (t[j+1] - t[j-1]) * t[j] for j in 2:length(t2YM[i,:])-1]
            global t_aux = t[2:end-1]
        else
            global t_aux = t[2:end-1]
            ixm = findmin(abs.(t_aux .- (w0-0.5)))[2]
            ixM = findmin(abs.(t_aux[1:end-1] .- (w0+0.5)))[2]
            tdt2YM[i,1:ixm-1] .= t2YM[i,2]
            tdt2YM[i,ixm:end] = [(t2YM[i,j+1] - t2YM[i,j-1]) / (t[j+1] - t[j-1]) * t[j] for j in ixm+1:length(t2YM[i,:])-1]
        end
    end

    return t2YM, tdt2YM, W_obs, t, t_aux
end

function concat_data!(data1::Vector{juobs.CData}, data2::Vector{juobs.CData})
    N = length(data1)
    if length(data1) != length(data2) 
        error("number of correlators do not match")
    end
    for k = 1:N
        data1[k].vcfg = vcat(data1[k].vcfg, data2[k].vcfg)
        data1[k].re_data = vcat(data1[k].re_data, data2[k].re_data)
        data1[k].im_data = vcat(data1[k].im_data, data2[k].im_data)
        idx = sortperm(data1[k].vcfg)
        data1[k].vcfg = data1[k].vcfg[idx]
        data1[k].re_data = data1[k].re_data[idx, :] 
        data1[k].im_data = data1[k].im_data[idx, :] 
    end
    return nothing
end

function concat_data!(data1::Vector{Vector{juobs.CData}}, data2::Vector{Vector{juobs.CData}})
    N = length(data1)
    if length(data1) != length(data2) 
        error("number of correlators do not match")
    end
    R = length(data1[1])
    if length(data1[1]) != length(data2[1])
        error("number of replicas do not match")
    end
    for k = 1:N
        for r = 1:R
            data1[k][r].vcfg = vcat(data1[k][r].vcfg, data2[k][r].vcfg)
            data1[k][r].re_data = vcat(data1[k][r].re_data, data2[k][r].re_data)
            data1[k][r].im_data = vcat(data1[k][r].im_data, data2[k][r].im_data)
            idx = sortperm(data1[k][r].vcfg)
            data1[k][r].vcfg = data1[k][r].vcfg[idx]
            data1[k][r].re_data = data1[k][r].re_data[idx, :] 
            data1[k][r].im_data = data1[k][r].im_data[idx, :] 

        end
    end
    return nothing
end


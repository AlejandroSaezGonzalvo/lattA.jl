using juobs

function get_corr_wil(path::String, ens::EnsInfo, g1::String, g2::String; rw=false, info=false, legacy=false, fs=false)
    path_rw = joinpath(path, ens.id, "rwf")
    path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    path = joinpath(path, ens.id, "wil")
    path = filter(x->occursin(".mesons.dat", x), readdir(path, join=true))

    rwf = read_ms1.(path_rw, v=ens.vrw)
    dat = read_mesons([path[i] for i in 1:length(path)], g1, g2, legacy=legacy)
    truncate_data!(dat,size.(rwf,2))

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

function get_YW(path::String, ens::EnsInfo; rw=false, npol::Int64=2, ws::ADerrors.wspace=ADerrors.wsg)

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
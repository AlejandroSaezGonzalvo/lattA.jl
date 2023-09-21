using juobs

function get_corr_wil(path::String, ens::EnsInfo, g1::String, g2::String; rw=false, info=false, legacy=false, fs=false)
    path_rw = joinpath(path, ens.id, "rwf")
    path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    path = joinpath(path, ens.id, "wil")
    path = filter(x->occursin(".mesons.dat", x), readdir(path, join=true))

    rwf = read_ms1.(path_rw)
    dat = read_mesons([path[i] for i in 1:length(path)], g1, g2, legacy=legacy)
    rw ? corr = [corr_obs(dat[i], L=ens.L, rw=rwf, info=info, flag_strange=fs) for i in 1:length(dat)] : corr = [corr_obs(dat[i], L=L[index], info=info, flag_strange=fs) for i in 1:length(dat)]

    if info == false
        return corr
    else
        pp = getindex.(corr,1)
        ppw = getindex.(corr,2)
        w = getindex.(corr,3)
        return pp, ppw, w[1]
    end
end
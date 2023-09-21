function get_corr(path::Vector{String}, ens::EnsInfo, g1::String, g2::String; rw=nothing, info=false, legacy=false, fs=false)
    path = joinpath.(path, ens.id)
    path = filter(x->occursin(".dat", x), readdir(p, join=true))

    dat = read_mesons([path[i] for i in 1:length(path)], g1, g2, legacy=legacy)
    corr = [corr_obs(dat[i], L=L[index], rw=rw, info=info, flag_strange=fs) for i in 1:length(dat)]

    if info == false
        return corr
    else
        pp = getindex.(corr,1)
        ppw = getindex.(corr,2)
        w = getindex.(corr,3)
        return pp, ppw, w[1]
    end
end
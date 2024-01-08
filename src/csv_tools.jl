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

    return juobs.Corr(obs, k, mu, gamma, y0)
end

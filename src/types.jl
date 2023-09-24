mutable struct EnsInfo
    id::String
    L::Int64
    beta::Float64
    ca::Float64
    dtr::Int64
    vrw::String
    function EnsInfo(ens_id::String, info::Vector{Any})
        id = ens_id
        L = info[1]
        beta = info[2]
        dtr = info[3]
        vrw = info[4]

        p0 = 9.2056
        p1 = -13.9847
        g2 = 6 ./ beta 
        ca = - 0.006033 .* g2 .*( 1 .+exp.(p0 .+ p1./g2)) 
        return new(id, L, beta, ca, dtr, vrw)
    end
end
    
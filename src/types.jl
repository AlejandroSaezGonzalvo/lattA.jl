mutable struct EnsInfo
    id::String
    L::Int64
    T::Int64
    beta::Float64
    ca::uwreal
    dtr::Int64
    vrw::Union{String, Vector{String}}
    cnfg::Vector{Int64}
    function EnsInfo(ens_id::String, info::Vector{Any})
        id = ens_id
        L = info[1]
        T = info[2]
        beta = info[3]
        dtr = info[4]
        vrw = info[5]
	    cnfg = info[6]

        p0 = 9.2056
        p1 = -13.9847
        g2 = 6 ./ beta 
        ca = - 0.006033 .* g2 .*( 1 .+exp.(p0 .+ p1./g2)) 
        err_ca = 0.08 * ca
        ca_uw = uwreal([ca, err_ca], string("ca-",ens_id))
        #return new(id, L, T, beta, ca, dtr, vrw, cnfg)
        return new(id, L, T, beta, ca_uw, dtr, vrw, cnfg)
    end
end
    

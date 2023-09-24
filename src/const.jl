ensemble = ["H101", "H102r001", "H102r002", "H105", "H105r005",
            "H400",
            "N200", "N202", "N203", "D200",
            "N300", "J303"]

const ens_db = Dict(
    #"ens_id"=>[L, beta, dtr, vrw]
    "H101"     => [32, 3.4, 2, "1.2"],
    "H102r001" => [32, 3.4, 2, "1.2"],
    "H102r002" => [32, 3.4, 2, "1.2"],
    "H105"     => [32, 3.4, 2, "1.2"],
    "H105r005" => [32, 3.4, 1, "1.2"],
    "H400"     => [32, 3.46, 1, "1.2"],
    "N200"     => [48, 3.55, 1, "1.2"],
    "N202"     => [48, 3.55, 2, "1.2"],
    "N203"     => [48, 3.55, 1, "1.4"],
    "D200"     => [64, 3.55, 2, "1.2"],
    "N300"     => [48, 3.70, 1, "1.2"],
    "N302"     => [],
    "J303"     => [64, 3.70, 2, "1.2"]
)

#PDG & FLAG21
const hc = 197.3269804 #MeV fm
const Mpi = uwreal([134.9768,0.0005],"mpi PDG") #MeV
const Mk = uwreal([497.611,0.013],"mk PDG") #Mev
const fpi = uwreal([130.56,0.02],"fpi exp") + uwreal([0.0,0.13],"fpi QED") + uwreal([0.0,0.02],"fpi Vud") #MeV FLAG21
const fk = uwreal([157.2,0.2],"fk exp") + uwreal([0.0,0.2],"fk QED") + uwreal([0.0,0.4],"fk Vud") #MeV FLAG21

const ens_obs = Dict(
    #"ens_id"=>[t0,mpi,mk,m12,m13,fpi,fk]
    "H101"     => [2.86, 0.18, 0.18, 0.0092, 0.0092, 0.083, 0.083]
)

wpm = Dict{String, Vector{Float64}}()
for e in ensemble wpm[e] = [-1.0, -1.0, 4.0, -1.0] end

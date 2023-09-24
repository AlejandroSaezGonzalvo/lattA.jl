const ensemble = Dict(
    1  => "H101",
    2  => "H102r001",
    3  => "H102r002",
    4  => "H105",
    5  => "H105r005",
    6  => "H400",
    7  => "N202",
    8  => "N203",
    9  => "N200",
    10 => "D200",
    11 => "N300",
    12 => "J303"
)

const ens_db = Dict(
    #"ens_id"=>[L, beta, dtr, vrw]
    "H101"     => [32, 3.4, 2, "1.2"],
    "H102r001" => [32, 3.4, 2, "1.2"],
    "H102r002" => [32, 3.4, 2, "1.2"],
    "H105"     => [32, 3.4, 2, "1.2"],
    "H105r005" => [32, 3.4, 1, "1.2"],
    "H400"     => [32, 3.46, 1, "1.2"],
    "N202"     => [48, 3.55, 2, "1.2"],
    "N203"     => [48, 3.55, 1, "1.4"],
    "N200"     => [48, 3.55, 1, "1.2"],
    "D200"     => [64, 3.55, 2, "1.2"],
    "N300"     => [48, 3.70, 1, "1.2"],
    "J303"     => [64, 3.70, 2, "1.2"]
)

wpm = Dict(
    #"ens_id"=>[L, beta, dtr, vrw]
    "H101"     => [-1.0, -1.0, 4.0, -1.0],
    "H102r001" => [-1.0, -1.0, 4.0, -1.0],
    "H102r002" => [-1.0, -1.0, 4.0, -1.0],
    "H105"     => [-1.0, -1.0, 4.0, -1.0],
    "H105r005" => [-1.0, -1.0, 4.0, -1.0],
    "H400"     => [-1.0, -1.0, 4.0, -1.0],
    "N202"     => [-1.0, -1.0, 4.0, -1.0],
    "N203"     => [-1.0, -1.0, 4.0, -1.0],
    "N200"     => [-1.0, -1.0, 4.0, -1.0],
    "D200"     => [-1.0, -1.0, 4.0, -1.0],
    "N300"     => [-1.0, -1.0, 4.0, -1.0],
    "J303"     => [-1.0, -1.0, 4.0, -1.0]
)

#PDG & FLAG21
const hc = 197.3269804 #MeV fm
const Mpi = uwreal([134.9768,0.0005],"mpi PDG") #MeV
const MK = uwreal([497.611,0.013],"mk PDG") #Mev
const Fpi = uwreal([130.56,0.02],"fpi exp") + uwreal([0.0,0.13],"fpi QED") + uwreal([0.0,0.02],"fpi Vud") #MeV FLAG21
const FK = uwreal([157.2,0.2],"fk exp") + uwreal([0.0,0.2],"fk QED") + uwreal([0.0,0.4],"fk Vud") #MeV FLAG21

const ens_obs = Dict(
    #"ens_id"=>[t0,mpi,mk,m12,m13,fpi,fk]
    "H101"     => [2.86, 0.18, 0.18, 0.0092, 0.0092, 0.083, 0.083],
    "H102r001" => [2.88, 0.15, 0.19, 0.0063, 0.0100, 0.078, 0.083],
    "H102r002" => [2.88, 0.15, 0.19, 0.0065, 0.0101, 0.080, 0.086],
    "D200"     => [5.18, 0.65, 0.15, 0.0015, 0.0094, 0.054, 0.064]
)


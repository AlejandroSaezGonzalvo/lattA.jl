using LaTeXStrings

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
    11 => "E250",
    12 => "N300",
    13 => "N302",
    14 => "J303",
    15 => "E300",
    16 => "J500",
    17 => "J501"
)

const ensemble_inv = Dict(
    "H101"     => 1,
    "H102r001" => 2,
    "H102r002" => 3,
    "H105"     => 4,
    "H105r005" => 5,
    "H400"     => 6,
    "N202"     => 7,
    "N203"     => 8,
    "N200"     => 9,
    "D200"     => 10,
    "E250"     => 11,
    "N300"     => 11,
    "N302"     => 13,
    "J303"     => 14,
    "E300"     => 15,
    "J500"     => 16,
    "J501"     => 17
)

const ens_db = Dict(
    #"ens_id"=>[L, T, beta, dtr, vrw, trunc_cnfg]
    "H101"     => [32, 96, 3.4, 2, "1.2", [1001,1009]],
    "H102r001" => [32, 96, 3.4, 2, "1.2", [997]],
    "H102r002" => [32, 96, 3.4, 2, "1.2", [1008]],
    "H105"     => [32, 96, 3.4, 2, "1.2", [947,1042]],
    "H105r005" => [32, 96, 3.4, 1, "1.2", [837]],
    "H400"     => [32, 96, 3.46, 1, "1.2", [505,540]],
    "N202"     => [48, 128, 3.55, 2, "1.2", [899]],
    "N203"     => [48, 128, 3.55, 1, "1.4", [756,787]],
    "N200"     => [48, 128, 3.55, 1, "1.2", [856,856]],
    "D200"     => [64, 128, 3.55, 2, "1.2", [2001]],
    "E250"     => [96, 192, 3.55, 1, ["1.6"], [1009]],
    "N300"     => [48, 128, 3.70, 1, "1.2", [1521]],
    "N302"     => [48, 128, 3.70, 1, "1.2", [2201]],
    "J303"     => [64, 192, 3.70, 2, "1.2", [1073]],
    "E300"     => [96, 192, 3.70, 1, ["1.4"], [1139]],
    "J500"     => [64, 192, 3.85, 2, ["1.2", "1.4", "2.0"], [789,655,431]],
    "J501"     => [64, 192, 3.85, 1, ["1.2", "1.4", "2.0"], [1635,1142,1150]]
)

const ens_cnfg = Dict(
    #"ens_id"=>[cnfg_rep_true]
    "H101"     => [1001,1009],
    "H102r001" => [997],
    "H102r002" => [1008],
    "H105"     => [1027,1042],
    "H105r005" => [837],
    "H400"     => [505,540],
    "N202"     => [899],
    "N203"     => [756, 787],
    "N200"     => [856, 856],
    "D200"     => [2001],
    "N300"     => [1521],
    "N302"     => [2201],
    "J303"     => [1073],
    "J501"     => [50]
)

const ens_kappa = Dict(
    #"ens_id"=>kappa_grid
    "H101"     => [0.137250, 0.137276, 0.137300],
    "H102r001" => [0.1372457, 0.1372954, 0.1373451],
    "H102r002" => [0.1372457, 0.1372954, 0.1373451],
    "H105"     => [0.1372613, 0.1372977, 0.1373340],
    "H105r005" => [0.1372613, 0.1372977, 0.1373340],
    "H400"     => [0.137272, 0.137300, 0.137330],
    "N202"     => [0.137287, 0.1372986, 0.137310],
    "N203"     => [0.137285, 0.1373057, 0.137325],
    "N200"     => [0.1372947, 0.137300, 0.1373249],
    "D200"     => [0.1373004, 0.1373132, 0.1373261],
    "E250"     => [0.1373072679, 0.137321],
    "N300"     => [0.137180, 0.137206, 0.137218],
    "N302"     => [0.137204, 0.137215],
    "E300"     => [0.1372059, 0.1372222],
    "J303"     => [0.1372036, 0.1372118, 0.1372199],
    "J500"     => [0.1369833003, 0.1370106997],
    "J501"     => [0.1369872999, 0.1370147001],
)

const ens_mul = Dict(
    #"ens_id"=>mul_grid
    "H101"     => [0.00645, 0.0067, 0.006945],
    "H102r001" => [0.004721, 0.004862, 0.005008],
    "H102r002" => [0.004721, 0.004862, 0.005008],
    "H105"     => [0.002800, 0.002920, 0.003040],
    "H105r005" => [0.002800, 0.002920, 0.003040],
    "H400"     => [0.005746, 0.006000, 0.006260],
    "N202"     => [0.005092, 0.005196, 0.005230],
    "N203"     => [0.003500, 0.003576, 0.003700],
    "N200"     => [0.002318, 0.002390, 0.002461],
    "D200"     => [0.001200, 0.001300],
    "E250"     => [0.000442, 0.0004784],
    "N300"     => [0.0037, 0.0041, 0.0044],
    "N302"     => [0.002628, 0.003212],
    "E300"     => [0.0006916, 0.0008284],
    "J303"     => [0.00135, 0.00165],
    "J500"     => [0.0029925, 0.0036575],
    "J501"     => [0.0019476, 0.0023804]
)

const ens_mus = Dict(
    #"ens_id"=>mus_grid
    "H101"     => [0.006450, 0.006700, 0.006945],
    "H102r001" => [0.010101, 0.010404, 0.010716],
    "H102r002" => [0.010101, 0.010404, 0.010716],
    "H105"     => [0.013200, 0.013700, 0.014200],
    "H105r005" => [0.013200, 0.013700, 0.014200],
    "H400"     => [0.005746, 0.006000, 0.006260],
    "N202"     => [0.005092, 0.005196, 0.005230],
    "N203"     => [0.008100, 0.008329, 0.008600],
    "N200"     => [0.010330, 0.010650, 0.010969],
    "D200"     => [0.012400, 0.013400],
    "E250"     => [0.013776, 0.014924],
    "N300"     => [0.0037, 0.0041, 0.0044],
    "N302"     => [0.006408, 0.007832],
    "E300"     => [0.01053, 0.01211],
    "J303"     => [0.0088, 0.0100],
    "J500"     => [0.0029925, 0.0036575],
    "J501"     => [0.005076, 0.006204]
)

const ens_up_fpik = Dict(
    #"ens_id"=>guess_up_fpik
    "H101"     => [-2530.6554321621106, 2.0171877597068244, 695.3460200403485, -47.762974774302045],
    "H102r001" => rand(5),
    "H102r002" => rand(5),
    "H105"     => [-3175.8070342116607, 2.0407450945085954, 0.600503353400361, 872.6430445629916, -59.943970541569044],
    "H105r005" => [-3175.8070342116607, 2.0407450945085954, 0.600503353400361, 872.6430445629916, -59.943970541569044],
    "H400"     => rand(5),
    "N202"     => [-22234.388160820763, 2.4770858675898504, 6106.142495506919, -419.22458075483473, -0.7772612475488999],
    "N203"     => [-12547.361232557383, 2.0023078576283284, 0.8099693106098493, 3446.0561021445487, -236.60768684371956],
    "N200"     => [-5095.750922147775, 2.237783076326602, 0.7458561664575111, 1399.9044528574202, -96.1435556034239],
    "D200"     => [-5869.087702588715, 2.372604551696622, 0.6491095784676267, 1612.7038244691962, -110.78251019157196],
    "E250"     => rand(5),
    "N300"     => [-99330.55686360868, 3.148674936123974, 27257.020477018694, -1869.878789249911, 1.660610917482043],
    "N302"     => rand(5),
    "J303"     => [-4907.782658533723, 2.6113901374563184, 0.8622480656404614, 1347.9014273915923, -92.54690971214653],
    "J500"     => rand(5),
    "J501"     => rand(5)
)

const ens_up_fpi = Dict(
    #"ens_id"=>guess_up_fpik
    "H101"     => rand(4),
    "H102r001" => rand(4),
    "H102r002" => rand(4),
    "H105"     => rand(4),
    "H105r005" => rand(4),
    "H400"     => rand(4),
    "N202"     => rand(4),
    "N203"     => [-6056.093837842525, 1.4662744627452966, 1663.3163101991515, -114.20734176958187],
    "N200"     => [-600.9070616956267, 1.545788321311486, 165.4174240331051, -11.383148177250593],
    "D200"     => [-16356.71950962457, 0.7030522823712142, 4492.653808265466, -308.49529516532726],
    "E250"     => rand(4),
    "N300"     => rand(4),
    "N302"     => rand(4),
    "J303"     => [-36488.99420061363, 1.3361050672576444, 10013.946389308961, -687.049848333641],
    "J500"     => rand(4),
    "J501"     => rand(4)
)

const ens_up_fk = Dict(
    #"ens_id"=>guess_up_fpik
    "H101"     => rand(5),
    "H102r001" => rand(5),
    "H102r002" => rand(5),
    "H105"     => rand(5),
    "H105r005" => rand(5),
    "H400"     => rand(5),
    "N202"     => rand(5),
    "N203"     => [-5262.492498397987, 0.5907616468097243, 0.5351661361257436, 1445.287897417873, -99.23242234315394],
    "N200"     => [-5105.184068785293, 0.7012124919983777, 0.4941984848363852, 1402.125796447806, -96.27174098418071],
    "D200"     => [-5869.290835968356, 1.2312353617579648, 0.42808942401543026, 1612.126763488127, -110.70048698560967],
    "E250"     => rand(5),
    "N300"     => rand(5),
    "N302"     => rand(5),
    "J303"     => [-4913.626695370824, 0.655423129636025, 0.44221644588819253, 1348.7216742914009, -92.55064727595374],
    "J500"     => rand(5),
    "J501"     => rand(5)
)

const ens_obs = Dict(
    #"ens_id"=>[t0,mpi,mk,m12,m13,fpi,fk] approx.
    "H101"     => [2.86, 0.18, 0.18, 0.0092, 0.0092, 0.083, 0.083],
    "H102r001" => [2.88, 0.15, 0.19, 0.0063, 0.0100, 0.078, 0.083],
    "H102r002" => [2.88, 0.15, 0.19, 0.0065, 0.0101, 0.080, 0.086],
    "D200"     => [5.18, 0.065, 0.15, 0.0015, 0.0094, 0.054, 0.064],
    "J501"     => [13.94, 0.066, 0.088, 0.0027, 0.0049, 0.035, 0.037]
)

const sym_bool = Dict(
    "J501" => true,
    "J500" => true,
    "E250" => false,
    "E300" => false
)

const db = Dict(
    "J501" => ["ts001_chunk1", "ts001_chunk3", "ts190_chunk1", "ts190_chunk2"],
    "J500" => ["ts001_chunk1", "ts001_chunk2", "ts190_chunk1", "ts190_chunk2"],
    "E250" => ["ts001", "ts097"],
    "E300" => ["ts001_chunk2", "ts001_chunk3", "ts190_chunk1"]#, "ts001_chunk1"]
)

const db_c = Dict(
    "J501" => ["ts001", "ts190"],
    "J500" => ["ts001_chunk1", "ts190_chunk2"],
    "E250" => ["ts001"],
    "E300" => ["ts001_chunk2"]#, "ts001_chunk1"]
)

const beta_ZA = Dict(
    #beta=>ZA
    3.40 => uwreal([0.75642,0.00072], "b=3.40"),
    3.46 => uwreal([0.76169,0.00093], "b=3.46"),
    3.55 => uwreal([0.76979,0.00043], "b=3.55"),
    3.70 => uwreal([0.78378,0.00047], "b=3.70"),
    3.85 => uwreal([0.79667,0.00047], "b=3.85")
)

const beta_ZP = Dict(
    #beta=>ZP
    3.40 => uwreal([0.35121,0.00056], "b=3.40"),
    3.46 => uwreal([0.34941,0.00044], "b=3.46"),
    3.55 => uwreal([0.34767,0.00055], "b=3.55"),
    3.70 => uwreal([0.34732,0.00063], "b=3.70"),
    3.85 => uwreal([0.35014,0.00073], "b=3.85")
)

const beta_bap = Dict(
    #beta=>ZP
    3.40 => uwreal([-0.324,0.017], "b=3.40"),
    3.46 => uwreal([-0.265,0.014], "b=3.46"),
    3.55 => uwreal([-0.196,0.014], "b=3.55"),
    3.70 => uwreal([-0.119,0.014], "b=3.70"),
    3.85 => uwreal([-0.073,0.012], "b=3.85")
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
    "E250"     => [-1.0, -1.0, 4.0, -1.0],
    "N300"     => [-1.0, -1.0, 4.0, -1.0],
    "N302"     => [-1.0, -1.0, 4.0, -1.0],
    "J303"     => [-1.0, -1.0, 4.0, -1.0],
    "E300"     => [-1.0, -1.0, 4.0, -1.0],
    "J500"     => [-1.0, -1.0, 4.0, -1.0],
    "J501"     => [-1.0, -1.0, 4.0, -1.0]
)

tm_mpi = Dict(
    "E300" => [[1], [40,50,60]],
    "J500" => [[1], [30,40,50,90]],
    "J501" => [[1], [20,30,40,50,60]]
)

tm_mk = Dict(
    "E300" => [[1], [90,100,110,120,140,145,150]],
    "J500" => [[1], [30,40,50,90]],
    "J501" => [[1], [40,50,60,70]]
)

tm_m12 = Dict(
    "E300" => [[1], [40,50,60]],
    "J500" => [[1], [30,50,60]],
    "J501" => [[1], [40,50,59]]
)

tm_m13 = Dict(
    "E300" => [[1], [80,90,100]],
    "J500" => [[1], [30,50,60]],
    "J501" => [[1], collect(60:68)]
)

tm_fpi = Dict(
    "E300" => [[1], 96 .- [10,15,20,25,40]],
    "J500" => [[1], [35,40,55]],
    "J501" => [[1], [40,50,60,70]]
)

tm_fk = Dict(
    "E300" => [[1], [50,60,80]],
    "J500" => [[1], [35,40,55]],
    "J501" => [[1], [30,40,50,60,70]]
)

tM_mpi = Dict(
    "E300" => [[180], [120,130,140]],
    "J500" => [[11], [100,120,130,140,160]],
    "J501" => [[180], [120,130,140,150,160]]
)

tM_mk = Dict(
    "E300" => [[11], [90,100,110,120,140,145,150]],
    "J500" => [[11], [100,120,130,140,160]],
    "J501" => [[180], [120,130,140,150,160]]
)

tM_m12 = Dict(
    "E300" => [[11], [120,130,140]],
    "J500" => [[11], [120,150,160]],
    "J501" => [[11], [120,130,140]]
)

tM_m13 = Dict(
    "E300" => [[11], [120,130,140]],
    "J500" => [[11], [120,150,160]],
    "J501" => [[11], collect(124:132)]
)

tM_fpi = Dict(
    "E300" => [[11], 96 .+ [10,15,20,25,40]],
    "J500" => [[11], [120,130,151,160]],
    "J501" => [[180], [120,130,140,150,160]]
)

tM_fk = Dict(
    "E300" => [[11], [110,124,135,150]],
    "J500" => [[11], [120,130,151,160]],
    "J501" => [[180], [120,130,140,150,160]]
)

function bAt(beta::Float64)
    return 1 + 0.0472 * (6 / beta)
end

#PDG & FLAG21
const hc = 197.3269804 #MeV fm
const Mpi = uwreal([134.9768,0.0005],"mpi PDG") #MeV
const MK = uwreal([497.611,0.013],"mk PDG") #Mev
const Fpi = uwreal([130.56,0.02],"fpi exp") + uwreal([0.0,0.13],"fpi QED") + uwreal([0.0,0.02],"fpi Vud") #MeV FLAG21
const FK = uwreal([157.2,0.2],"fk exp") + uwreal([0.0,0.2],"fk QED") + uwreal([0.0,0.4],"fk Vud") #MeV FLAG21

#t0 model av model tags
const mods = [L"$[\chi SU(3)][a^2][-]$", 
        L"$[\chi SU(3)][a^2][\beta>3.40]$",
        L"$[\chi SU(3)][a^2][\beta>3.46]$", 
        L"$[\chi SU(3)][a^2][m_{\pi}<420\;MeV]$",
        L"$[\chi SU(3)][a^2][\beta>3.40\;&\;m_{\pi}<420\;MeV]$",
        L"$[\chi SU(3)][a^2][m_{\pi}<350\;MeV]$", 
        L"$[\chi SU(3)][a^2][m_{\pi}L>3.9]$",
        L"$[\chi SU(3)][a^2][m_{\pi}L>4.1]$", 
        L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][-]$", 
        L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][\beta>3.40]$",
        L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][\beta>3.46]$", 
        L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}<420\;MeV]$",
        L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][\beta>3.40\;&\;m_{\pi}<420\;MeV]$",
        L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}<350\;MeV]$",
        L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}L>3.9]$",
        L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}L>4.1]$",
        L"$[\chi SU(3)][a^2+a^2\phi_2][-]$", 
        L"$[\chi SU(3)][a^2+a^2\phi_2][\beta>3.40]$", 
        L"$[\chi SU(3)][a^2+a^2\phi_2][\beta>3.46]$", 
        L"$[\chi SU(3)][a^2+a^2\phi_2][m_{\pi}<420\;MeV]$",
        L"$[\chi SU(3)][a^2+a^2\phi_2][\beta>3.40\;&\;m_{\pi}<420\;MeV]$",
        L"$[\chi SU(3)][a^2+a^2\phi_2][m_{\pi}<350\;MeV]$",
        L"$[\chi SU(3)][a^2+a^2\phi_2][m_{\pi}L>3.9]$",
        L"$[\chi SU(3)][a^2+a^2\phi_2][m_{\pi}L>4.1]$", 
        L"$[Tay][a^2][-]$", 
        L"$[Tay][a^2][\beta>3.40]$", 
        L"$[Tay][a^2][\beta>3.46]$", 
        L"$[Tay][a^2][m_{\pi}<420\;MeV]$", 
        L"$[Tay][a^2][\beta>3.40\;&\;m_{\pi}<420\;MeV]$",
        L"$[Tay][a^2][m_{\pi}<350\;MeV]$", 
        L"$[Tay][a^2][m_{\pi}L>3.9]$",
        L"$[Tay][a^2][m_{\pi}L>4.1]$", 
        L"$[Tay4][a^2][-]$", 
        L"$[Tay4][a^2][\beta>3.40]$", 
        L"$[Tay4][a^2][\beta>3.46]$", 
        L"$[Tay4][a^2][m_{\pi}<420\;MeV]$",
        L"$[Tay4][a^2][\beta>3.40\;&\;m_{\pi}<420\;MeV]$", 
        L"$[Tay4][a^2][m_{\pi}<350\;MeV]$", 
        L"$[Tay4][a^2][m_{\pi}L>3.9]$",
        L"$[Tay4][a^2][m_{\pi}L>4.1]$", 
        L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][-]$", 
        L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][\beta>3.40]$", 
        L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][\beta>3.46]$", 
        L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}<420\;MeV]$",
        L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][\beta>3.40\;&\;m_{\pi}<420\;MeV]$",
        L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}<350\;MeV]$",
        L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}L>3.9]$",
        L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}L>4.1]$", 
        L"$[Tay][a^2+a^2\phi_2][-]$", 
        L"$[Tay][a^2+a^2\phi_2][\beta>3.40]$", 
        L"$[Tay][a^2+a^2\phi_2][\beta>3.46]$", 
        L"$[Tay][a^2+a^2\phi_2][m_{\pi}<420\;MeV]$",
        L"$[Tay][a^2+a^2\phi_2][\beta>3.40\;&\;m_{\pi}<420\;MeV]$",
        L"$[Tay][a^2+a^2\phi_2][m_{\pi}<350\;MeV]$",
        L"$[Tay][a^2+a^2\phi_2][m_{\pi}L>3.9]$",
        L"$[Tay][a^2+a^2\phi_2][m_{\pi}L>4.1]$",
        L"$[\chi SU(2)][a^2][-]$", 
        L"$[\chi SU(2)][a^2][\beta>3.40]$",
        L"$[\chi SU(2)][a^2][\beta>3.46]$", 
        L"$[\chi SU(2)][a^2][m_{\pi}<420\;MeV]$",
        L"$[\chi SU(2)][a^2][\beta>3.40\;&\;m_{\pi}<420\;MeV]$",
        L"$[\chi SU(2)][a^2][m_{\pi}<350\;MeV]$", 
        L"$[\chi SU(2)][a^2][m_{\pi}L>3.9]$",
        L"$[\chi SU(2)][a^2][m_{\pi}L>4.1]$"]

const mods_346 = [L"$[\chi SU(3)][a^2][-]$", 
        #L"$[\chi SU(3)][a^2][\beta>3.40]$",
        L"$[\chi SU(3)][a^2][\beta>3.46]$", 
        L"$[\chi SU(3)][a^2][m_{\pi}<420\;MeV]$",
        #L"$[\chi SU(3)][a^2][\beta>3.40\;&\;m_{\pi}<420\;MeV]$",
        L"$[\chi SU(3)][a^2][m_{\pi}<350\;MeV]$", 
        L"$[\chi SU(3)][a^2][m_{\pi}L>3.9]$",
        L"$[\chi SU(3)][a^2][m_{\pi}L>4.1]$", 
        L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][-]$", 
        #L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][\beta>3.40]$",
        L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][\beta>3.46]$", 
        L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}<420\;MeV]$",
        #L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][\beta>3.40\;&\;m_{\pi}<420\;MeV]$",
        L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}<350\;MeV]$",
        L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}L>3.9]$",
        L"$[\chi SU(3)][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}L>4.1]$",
        L"$[\chi SU(3)][a^2+a^2\phi_2][-]$", 
        #L"$[\chi SU(3)][a^2+a^2\phi_2][\beta>3.40]$", 
        L"$[\chi SU(3)][a^2+a^2\phi_2][\beta>3.46]$", 
        L"$[\chi SU(3)][a^2+a^2\phi_2][m_{\pi}<420\;MeV]$",
        #L"$[\chi SU(3)][a^2+a^2\phi_2][\beta>3.40\;&\;m_{\pi}<420\;MeV]$",
        L"$[\chi SU(3)][a^2+a^2\phi_2][m_{\pi}<350\;MeV]$",
        L"$[\chi SU(3)][a^2+a^2\phi_2][m_{\pi}L>3.9]$",
        L"$[\chi SU(3)][a^2+a^2\phi_2][m_{\pi}L>4.1]$", 
        L"$[Tay][a^2][-]$", 
        #L"$[Tay][a^2][\beta>3.40]$", 
        L"$[Tay][a^2][\beta>3.46]$", 
        L"$[Tay][a^2][m_{\pi}<420\;MeV]$", 
        #L"$[Tay][a^2][\beta>3.40\;&\;m_{\pi}<420\;MeV]$",
        L"$[Tay][a^2][m_{\pi}<350\;MeV]$", 
        L"$[Tay][a^2][m_{\pi}L>3.9]$",
        L"$[Tay][a^2][m_{\pi}L>4.1]$", 
        L"$[Tay4][a^2][-]$", 
        #L"$[Tay4][a^2][\beta>3.40]$", 
        L"$[Tay4][a^2][\beta>3.46]$", 
        L"$[Tay4][a^2][m_{\pi}<420\;MeV]$",
        #L"$[Tay4][a^2][\beta>3.40\;&\;m_{\pi}<420\;MeV]$", 
        L"$[Tay4][a^2][m_{\pi}<350\;MeV]$", 
        L"$[Tay4][a^2][m_{\pi}L>3.9]$",
        L"$[Tay4][a^2][m_{\pi}L>4.1]$", 
        L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][-]$", 
        #L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][\beta>3.40]$", 
        L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][\beta>3.46]$", 
        L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}<420\;MeV]$",
        #L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][\beta>3.40\;&\;m_{\pi}<420\;MeV]$",
        L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}<350\;MeV]$",
        L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}L>3.9]$",
        L"$[Tay][a^2\alpha_s^{\hat{\Gamma}}][m_{\pi}L>4.1]$", 
        L"$[Tay][a^2+a^2\phi_2][-]$", 
        #L"$[Tay][a^2+a^2\phi_2][\beta>3.40]$", 
        L"$[Tay][a^2+a^2\phi_2][\beta>3.46]$", 
        L"$[Tay][a^2+a^2\phi_2][m_{\pi}<420\;MeV]$",
        #L"$[Tay][a^2+a^2\phi_2][\beta>3.40\;&\;m_{\pi}<420\;MeV]$",
        L"$[Tay][a^2+a^2\phi_2][m_{\pi}<350\;MeV]$",
        L"$[Tay][a^2+a^2\phi_2][m_{\pi}L>3.9]$",
        L"$[Tay][a^2+a^2\phi_2][m_{\pi}L>4.1]$",
        L"$[\chi SU(2)][a^2][-]$", 
        #L"$[\chi SU(2)][a^2][\beta>3.40]$",
        L"$[\chi SU(2)][a^2][\beta>3.46]$", 
        L"$[\chi SU(2)][a^2][m_{\pi}<420\;MeV]$",
        #L"$[\chi SU(2)][a^2][\beta>3.40\;&\;m_{\pi}<420\;MeV]$",
        L"$[\chi SU(2)][a^2][m_{\pi}<350\;MeV]$", 
        L"$[\chi SU(2)][a^2][m_{\pi}L>3.9]$",
        L"$[\chi SU(2)][a^2][m_{\pi}L>4.1]$"]

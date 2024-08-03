#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO, PyPlot, LinearAlgebra, SpecialFunctions
using ADerrors: err

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
#include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/plot.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/IO_BDIO.jl")

#ens_av = ["A654", "H101", "H102", "H105", "N101", "C101", "C102", "D150", "B450", "S400", "N451", "D450", "D451", "D452", "N202", "N203", "N200", "D251", "D200", "D201", "E250", "N300", "N302", "J303", "J304", "E300", "J500", "J501"]

beta = [3.40,3.40,3.40,3.40,3.40,3.40,3.46,3.46,3.46,3.46,3.46,3.46,3.55,3.55,3.55,3.55,3.55,3.55,3.55,3.70,3.70,3.70,3.70,3.85,3.85]

#ens_Mz = ["H101", "H102", "N101", "C101", "B450", "S400", "N451", "D450", "N202", "N203", "N200", "D251", "D200", "E250", "N300", "J303", "J304", "E300"]
#ens_av_m = ens_Mz

#=
fb = BDIO_open("/home/asaez/cls_ens/results/w0_ratio.bdio", "r")
BDIO_seek!(fb)
w0_ph_ini = read_uwreal(fb)
=#

sqrt_t0_ph_ini = uwreal([0.1449,0.0008], "t0_ph ini")
t0_ph_ini = sqrt_t0_ph_ini ^ 2
w0_ph_ini = uwreal([0.1725,0.0005], "w0_ph ini")

plat = "cls"
path = "/home/asaez/cls_ens/codes/lattA.jl/plots/"
Simon = true
global ph_point = "t0"

wpm = Dict(
    #"ens_id"=>[L, beta, dtr, vrw]
    "A653"     => [-1.0, -1.0, 0.1, -1.0],
    "A654"     => [-1.0, -1.0, 0.1, -1.0],
    "H101"     => [-1.0, -1.0, 0.1, -1.0],
    "H102r001" => [-1.0, -1.0, 0.1, -1.0],
    "H102r002" => [-1.0, -1.0, 0.1, -1.0],
    "H105"     => [-1.0, -1.0, 0.1, -1.0],
    "H105r005" => [-1.0, -1.0, 0.1, -1.0],
    "N101"     => [-1.0, -1.0, 0.1, -1.0],
    "C101"     => [-1.0, -1.0, 0.1, -1.0],
    "C102"     => [-1.0, -1.0, 0.1, -1.0],
    "D150"     => [-1.0, -1.0, 0.1, -1.0],
    "B450"     => [-1.0, -1.0, 0.1, -1.0],
    "S400"     => [-1.0, -1.0, 0.1, -1.0],
    "N452"     => [-1.0, -1.0, 0.1, -1.0],
    "N451"     => [-1.0, -1.0, 0.1, -1.0],
    "H400"     => [-1.0, -1.0, 0.1, -1.0],
    "D450"     => [-1.0, -1.0, 0.1, -1.0],
    "D451"     => [-1.0, -1.0, 0.1, -1.0],
    "D452"     => [-1.0, -1.0, 0.1, -1.0],
    "N202"     => [-1.0, -1.0, 0.1, -1.0],
    "N203"     => [-1.0, -1.0, 0.1, -1.0],
    "N200"     => [-1.0, -1.0, 0.1, -1.0],
    "D200"     => [-1.0, -1.0, 0.1, -1.0],
    "D201"     => [-1.0, -1.0, 0.1, -1.0],
    "D251"     => [-1.0, -1.0, 0.1, -1.0],
    "E250"     => [-1.0, -1.0, 0.1, -1.0],
    "N300"     => [-1.0, -1.0, 0.1, -1.0],
    "N302"     => [-1.0, -1.0, 0.1, -1.0],
    "J303"     => [-1.0, -1.0, 0.1, -1.0],
    "J304"     => [-1.0, -1.0, 0.1, -1.0],
    "E300"     => [-1.0, -1.0, 0.1, -1.0],
    "J500"     => [-1.0, -1.0, 0.1, -1.0],
    "J501"     => [-1.0, -1.0, 0.1, -1.0]
)

LambdaQCD = 340 ## in Mev
Lam = LambdaQCD / hc
t0_old = t0_ph_ini

Gamma_0 = 0.0
Gamma_1 = -0.111
Gamma_2 = 0.247
Gamma_3 = 0.519
Gamma_4 = 0.668
Gamma_5 = 0.760

#=========== 2203.08676 ===============#

    w0_Mz = [
        uwreal([3.694,0.009], "H101 Mainz"),
        uwreal([3.759,0.011], "H102 Mainz"),
        uwreal([3.779,0.014], "H105 Mainz"),
        uwreal([3.789,0.005], "N101 Mainz"),
        uwreal([3.836,0.006], "C101 Mainz"),
        uwreal([4.845,0.018], "B450 Mainz"),
        uwreal([4.895,0.015], "S400 Mainz"),
        uwreal([4.896,0.005], "N451 Mainz"),
        uwreal([4.933,0.006], "D450 Mainz"),
        uwreal([7.020,0.029], "N202 Mainz"),
        uwreal([6.982,0.015], "N203 Mainz"),
        uwreal([7.040,0.014], "N200 Mainz"),
        uwreal([7.099,0.007], "D200 Mainz"),
        uwreal([7.176,0.006], "E250 Mainz"),
        uwreal([11.821,0.048], "N300 Mainz"),
        uwreal([11.784,0.049], "N302 Mainz"),
        uwreal([12.101,0.031], "J303 Mainz"),
        uwreal([12.163,0.015], "E300 Mainz")
    ]

    t0_Mz = [
        uwreal([2.847,0.005], "H101 Mainz"),
        uwreal([2.881,0.006], "H102 Mainz"),
        uwreal([2.889,0.007], "H105 Mainz"),
        uwreal([2.890,0.002], "N101 Mainz"),
        uwreal([2.913,0.003], "C101 Mainz"),
        uwreal([3.663,0.008], "B450 Mainz"),
        uwreal([3.686,0.007], "S400 Mainz"),
        uwreal([3.682,0.002], "N451 Mainz"),
        uwreal([3.696,0.003], "D450 Mainz"),
        uwreal([5.165,0.012], "N202 Mainz"),
        uwreal([5.146,0.006], "N203 Mainz"),
        uwreal([5.164,0.006], "N200 Mainz"),
        uwreal([5.179,0.003], "D200 Mainz"),
        uwreal([5.203,0.002], "E250 Mainz"),
        uwreal([8.544,0.019], "N300 Mainz"),
        uwreal([8.526,0.019], "N302 Mainz"),
        uwreal([8.621,0.010], "J303 Mainz"),
        uwreal([8.622,0.006], "E300 Mainz")
    ]

    mpi_Mz = [
        uwreal([0.18356,0.00048], "H101 Mainz"),
        uwreal([0.15457,0.00055], "H102 Mainz"),
        uwreal([0.12353,0.000128], "H105 Mainz"),
        uwreal([0.12236,0.00047], "N101 Mainz"),
        uwreal([0.09616,0.00065], "C101 Mainz"),
        uwreal([0.16108,0.00043], "B450 Mainz"),
        uwreal([0.13592,0.00043], "S400 Mainz"),
        uwreal([0.11089,0.00029], "N451 Mainz"),
        uwreal([0.08362,0.00039], "D450 Mainz"),
        uwreal([0.13423,0.00030], "N202 Mainz"),
        uwreal([0.11266,0.00023], "N203 Mainz"),
        uwreal([0.09234,0.00028], "N200 Mainz"),
        uwreal([0.06515,0.00028], "D200 Mainz"),
        uwreal([0.04217,0.00028], "E250 Mainz"),
        uwreal([0.10618,0.00024], "N300 Mainz"),
        uwreal([0.08725,0.00024], "N302 Mainz"),
        uwreal([0.06481,0.00023], "J303 Mainz"),
        uwreal([0.04367,0.00016], "E300 Mainz")
    ]

    mk_Mz = [
        uwreal([0.18356,0.00048], "H101 Mainz"),
        uwreal([0.19164,0.00048], "H102 Mainz"),
        uwreal([0.20251,0.00087], "H105 Mainz"),
        uwreal([0.20189,0.00028], "N101 Mainz"),
        uwreal([0.20578,0.00032], "C101 Mainz"),
        uwreal([0.16108,0.00043], "B450 Mainz"),
        uwreal([0.17056,0.00038], "S400 Mainz"),
        uwreal([0.17827,0.00018], "N451 Mainz"),
        uwreal([0.18393,0.00018], "D450 Mainz"),
        uwreal([0.13421,0.00029], "N202 Mainz"),
        uwreal([0.14413,0.00019], "N203 Mainz"),
        uwreal([0.15076,0.00021], "N200 Mainz"),
        uwreal([0.15615,0.00016], "D200 Mainz"),
        uwreal([0.15924,0.00008], "E250 Mainz"),
        uwreal([0.10618,0.00024], "N300 Mainz"),
        uwreal([0.11373,0.00032], "N302 Mainz"),
        uwreal([0.11964,0.00020], "J303 Mainz"),
        uwreal([0.12372,0.00013], "E300 Mainz")
    ]

    mO = [
        uwreal([0.6169,0.0012], "H101 Darmstadt"),
        uwreal([0.6422,0.0013], "H102 Darmstadt"),
        uwreal([0.6631,0.0013], "N101 Darmstadt"),
        uwreal([0.6794,0.0014], "C101 Darmstadt"),
        uwreal([0.5549,0.0028], "B450 Darmstadt"),
        uwreal([0.5788,0.0016], "S400 Darmstadt"),
        uwreal([0.5946,0.0007], "N451 Darmstadt"),
        uwreal([0.6111,0.0007], "D450 Darmstadt"),
        uwreal([0.6164,0.0008], "D452 Darmstadt"),
        uwreal([0.4649,0.0021], "N202 Darmstadt"),
        uwreal([0.4882,0.0011], "N203 Darmstadt"),
        uwreal([0.5068,0.0014], "N200 Darmstadt"),
        uwreal([0.5218,0.0006], "D200 Darmstadt"),
        uwreal([0.5285,0.0003], "E250 Darmstadt"),
        uwreal([0.3676,0.0023], "N300 Darmstadt"),
        uwreal([0.4023,0.0014], "J303 Darmstadt"),
        uwreal([0.4138,0.0006], "E300 Darmstadt"),
        uwreal([0.2851,0.0018], "J500 Darmstadt"),
        uwreal([0.2996,0.0008], "J501 Darmstadt")
    ]

    uwerr.(mO)
    uwerr.(w0_Mz)
    uwerr.(t0_Mz)
    uwerr.(mpi_Mz)
    uwerr.(mk_Mz)

    function read_data(obs::String="r")

        if obs == "omega"
            global ens_av_m = ["H101", "H102", "N101", "C101", "B450", "S400", "N451", "D450", "D452", "N202", "N203", "N200", "D200", "E250", "N300", "J303", "E300", "J500", "J501"]
        elseif obs == "r"
            global ens_av_m = ["H101", "H102", "N101", "C101", "C102", "D150", "B450", "S400", "N451", "D450", "D451", "D452", "N202", "N203", "N200", "D251", "D200", "D201", "E250", "N300", "J303", "J304", "E300", "J500", "J501"]
        end

        global ens_nosym = ["H102", "N101", "C101", "C102", "D150", "S400", "N451", "D450", "D451", "D452", "N203", "N200", "D251", "D200", "D201", "E250", "J303", "J304", "E300", "J501"]
        global ens_334 = findall(x -> x in ["A653", "A654"], ens_av_m)
        global ens_340 = findall(x -> x in ["H101", "H102", "H105", "N101", "C101", "C102", "D150"], ens_av_m)
        global ens_346 = findall(x -> x in ["B450", "S400", "N451", "D450", "D451", "D452"], ens_av_m)
        global ens_355 = findall(x -> x in ["N202", "N203", "N200", "D251", "D200", "D201", "E250"], ens_av_m)
        global ens_370 = findall(x -> x in ["N300", "N302", "J303", "J304", "E300"], ens_av_m)
        global ens_385 = findall(x -> x in ["J500", "J501"], ens_av_m)
        global ind_nosym = findall(x -> x in ens_nosym, ens_av_m)

        global w0 = Array{uwreal,1}()
        global t0 = Array{uwreal,1}()
        global w0_2 = Array{uwreal,1}()
        global t0_2 = Array{uwreal,1}()
        global w0_3 = Array{uwreal,1}()
        global t0_3 = Array{uwreal,1}()
        for id in ens_av_m
            obs = Array{uwreal,1}()
            obs_2 = Array{uwreal,1}()
            obs_3 = Array{uwreal,1}()
            if plat == "cls" && Simon == false
                fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted_w0/cls_w0_t0_", id, "_obs_wil_un.bdio"), "r") 
                BDIO_seek!(fb); push!(obs, read_uwreal(fb)); BDIO_seek!(fb, 2); push!(obs, read_uwreal(fb))
                BDIO_close!(fb)
                push!(w0,obs[1])
                push!(t0,obs[2])
            elseif Simon == true
                fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted_w0/cls_w0_t0_tau0.0_", id, "_obs_wil_un.bdio"), "r") 
                BDIO_seek!(fb); push!(obs, read_uwreal(fb)); BDIO_seek!(fb, 2); push!(obs, read_uwreal(fb))
                BDIO_close!(fb)
                push!(w0,obs[1])
                push!(t0,obs[2])
                fb_2 = BDIO_open(string("/home/asaez/cls_ens/results/unshifted_w0/cls_w0_t0_tau0.2_", id, "_obs_wil_un.bdio"), "r") 
                fb_3 = BDIO_open(string("/home/asaez/cls_ens/results/unshifted_w0/cls_w0_t0_tau-0.2_", id, "_obs_wil_un.bdio"), "r") 
                BDIO_seek!(fb_2); push!(obs_2, read_uwreal(fb_2)); BDIO_seek!(fb_2, 2); push!(obs_2, read_uwreal(fb_2))
                BDIO_close!(fb_2)
                push!(w0_2,obs_2[1])
                push!(t0_2,obs_2[2])
                BDIO_seek!(fb_3); push!(obs_3, read_uwreal(fb_3)); BDIO_seek!(fb_3, 2); push!(obs_3, read_uwreal(fb_3))
                BDIO_close!(fb_3)
                push!(w0_3,obs_3[1])
                push!(t0_3,obs_3[2])
            else 
                fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted_w0/new_w0_t0_", id, "_obs_wil_un.bdio"), "r") 
                BDIO_seek!(fb); push!(obs, read_uwreal(fb)); BDIO_seek!(fb, 2); push!(obs, read_uwreal(fb))
                BDIO_close!(fb)
                push!(w0,obs[1])
                push!(t0,obs[2])
            end
        end

        if Simon == true
            t0_p = similar(t0)
            w0_p = similar(w0)
            for i in 1:length(t0)
                t0_p[i] = t0[i] * w0_3[i] ^ 2 / w0_2[i] ^ 2
                w0_p[i] = w0[i] * t0_2[i] ^ 2 / t0[i] ^ 2
            end
            t0_or = deepcopy(t0)
            t0 = deepcopy(t0_p)
            w0 = deepcopy(w0_p)
        end

        global mpi = Array{uwreal,1}()
        global mk = Array{uwreal,1}()
        global t0_Ales = Array{uwreal,1}()
        for id in ens_av_m
            path_m = string("/home/asaez/cls_ens/results/unshifted_w0/", id, "_spectrum.bdio")
            aux = read_BDIO(path_m, "spectrum", "mpi")[1]
            push!(mpi, aux)
            aux = read_BDIO(path_m, "spectrum", "mk")[1]
            push!(mk, aux)
            aux = read_BDIO(path_m, "spectrum", "t0")[1]
            push!(t0_Ales, aux)
        end

        global phi2_Mz = 8 * t0_Mz .* mpi_Mz .^ 2
        global phi4_Mz = 8 * t0_Mz .* (mk_Mz .^ 2 .+ 0.5 * mpi_Mz .^ 2)
        global dphi = 8 * t0 .* (mk .^ 2 .- mpi .^ 2)
        global phi4star = 1.11
        global mO_ph = uwreal([1672.43,0.32], "mO exp") ## PDG18
        
        if ph_point == "t0"
            global phi2 = 8 * t0_or .* mpi .^ 2
            global phi4 = 8 * t0_or .* (mk .^ 2 .+ 0.5 * mpi .^ 2)
            global phi2_ph = 8 * t0_ph_ini * Mpi ^ 2 / hc ^ 2
            global phi4_ph = 8 * t0_ph_ini * (MK ^ 2 + 0.5 * Mpi ^ 2) / hc ^ 2
        elseif ph_point == "w0"
            global phi2 = 8 * w0 .* mpi .^ 2
            global phi4 = 8 * w0 .* (mk .^ 2 .+ 0.5 * mpi .^ 2)
            global Mpi = uwreal([134.9768,0.0005], "mpi input")
            global Mss = uwreal([689.89,0.49], "mss input")
            global MK = sqrt(0.5 * (Mss ^ 2 + Mpi ^ 2))
            global phi2_ph = 8 * w0_ph_ini ^ 2 * Mpi ^ 2 / hc ^ 2
            global phi4_ph = 8 * w0_ph_ini ^ 2 * (MK ^ 2 + 0.5 * Mpi ^ 2) / hc ^ 2
        end

        global wphi2 = 8 * w0 .* mpi .^ 2; uwerr.(wphi2)
        global wphiK = 8 * w0 .* mk .^ 2; uwerr.(wphiK)
        global wphi4 = wphiK .+ 0.5 * wphi2; uwerr.(wphi4)
        if obs == "omega" global wphiO = 8 * w0 .* mO .^ 2; uwerr.(mO) end
        global wphi2_ph = 8 * w0_ph_ini ^ 2 * Mpi ^ 2 / hc ^ 2
        global wphiK_ph = 8 * w0_ph_ini ^ 2 * MK ^ 2 / hc ^ 2
        global wphi4_ph = 8 * w0_ph_ini ^ 2 * (MK ^ 2 + 0.5 * Mpi ^ 2) / hc ^ 2

        global ratio1_ph = mO_ph ^ 2 / Mpi^2 / 8
        global ratio2_ph = mO_ph ^ 2 / (MK^2 + 0.5 * Mpi^2) / 8

        uwerr.(phi2)
        uwerr.(phi4)
        uwerr.(phi2_Mz)
        uwerr.(phi4_Mz)
        uwerr.(dphi)
        uwerr.(mpi)
        uwerr.(mk)
        uwerr.(w0)
        uwerr.(t0)
        uwerr.(t0_Ales)
    end

#===========================================#

#============ ratio ===================#

    read_data("r")

    #list = [7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25] ## no b=3.40
    #list = [2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,21,22,23,25] ## no sym
    #list = [8,9,10,11,12,14,15,16,17,18,19,21,22,23,25] ## no b=3.40 and no sym
    #list = [13,14,15,16,17,18,19,20,21,22,23,24,25] ## no b=3.46
    #list = [7,8,9,10,11,12,13,14,15,16,17,18,20,22,23,24,25] ## mpiL>4.1
    #list = [1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22,23,24,25] ## ??
    
    #list = [1,5,6,7,8,9,10,12,14,15,16,17,19,21,23,25]
    #list = [2,3,4,5,6,14,15,16,17,18,19,21,22,23,25]
    #list = [7,8,9,10,12,13,14,15,16,17,19,20,21,23,25]
    #list = [1,2,3,4,6,13,14,15,16,17,19,20,21,23,25]
    #list = [8,9,10,12,14,15,16,17,19,21,23,25]
    #list = ind_nosym

    y_Mz = sqrt.(t0_Mz) ./ sqrt.(w0_Mz)
    uwerr.(y_Mz)
    x_Mz = 1 ./ t0_Mz
    uwerr.(x_Mz)

    y = sqrt.(t0) ./ sqrt.(w0)
    uwerr.(y)
    x = [1 ./ t0 phi2 phi4]
    uwerr.(x)

    TIC_r = Array{Float64,1}()
    pval_r = Array{Float64,1}()
    dof_r = Array{Float64,1}()
    r_ph = Array{uwreal,1}()
    w0_ph = Array{uwreal,1}()

    ## funs
        function model_ch_a2as0(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as0_a3(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as0(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as0_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] + p[8] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as0_a3_a2phi2(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as0_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as0_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] + p[8] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as0_a3_a2phi4(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as0_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as1(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as1_a3(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as1(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as1_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] + p[8] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as1_a3_a2phi2(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as1_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as1_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] + p[8] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as1_a3_a2phi4(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as1_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as2_a3(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as2_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] + p[8] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as2_a3_a2phi2(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as2_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as2_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] + p[8] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as2_a3_a2phi4(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as2_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as3(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as3_a3(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as3(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as3_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] + p[8] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as3_a3_a2phi2(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as3_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as3_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] + p[8] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as3_a3_a2phi4(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as3_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as4_a3(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as4_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] + p[8] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as4_a3_a2phi2(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as4_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as4_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] + p[8] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as4_a3_a2phi4(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as4_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as5(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as5_a3(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as5(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as5_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] + p[8] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as5_a3_a2phi2(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as5_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as5_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] + p[8] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2as5_a3_a2phi4(x,p)
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2as5_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end





















        function model_Tay_a2as0(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as0_a3(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as0(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as0_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] + p[8] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as0_a3_a2phi2(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as0_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as0_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] + p[8] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as0_a3_a2phi4(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as0_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as1(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as1_a3(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as1(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as1_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] + p[8] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as1_a3_a2phi2(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as1_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as1_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] + p[8] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as1_a3_a2phi4(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as1_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as2_a3(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as2_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] + p[8] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as2_a3_a2phi2(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as2_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as2_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] + p[8] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as2_a3_a2phi4(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as2_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as3(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as3_a3(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as3(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as3_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] + p[8] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as3_a3_a2phi2(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as3_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as3_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] + p[8] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as3_a3_a2phi4(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as3_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as4_a3(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as4_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] + p[8] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as4_a3_a2phi2(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as4_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as4_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] + p[8] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as4_a3_a2phi4(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as4_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as5(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as5_a3(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as5(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as5_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] + p[8] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as5_a3_a2phi2(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as5_a2phi2(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] + p[9] * x[i,2] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as5_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] + p[8] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2as5_a3_a2phi4(x,p)
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else 
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2as5_a2phi4(x,p) 
            if x[1,1] == 0.0
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] + p[9] * x[i,3] * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end








        function model_ch_a2_a2as1_a3(x,p) 
            if x[1,1] == 0.0
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) for i in 1:div(length(x[:,1]),4)]
            else
                f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[9] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            end
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2(x,p)
            f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a3(x,p)
            f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a3_a4(x,p)
            f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] * x[i,1] ^ 0.5 + p[9] * x[i,1]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2mpi(x,p)
            f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2mpi_a2mk(x,p)
            f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] * x[i,2] + p[9] * x[i,3]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2mpi_a3(x,p)
            f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] * x[i,2] + p[9] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2mpi_a2mk_a3(x,p)
            f = [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] * x[i,2] + p[9] * x[i,3] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2(x,p)
            f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a3(x,p)
            f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a3_a4(x,p)
            f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] * x[i,1] ^ 0.5 + p[9] * x[i,1]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2mpi(x,p)
            f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+8] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2mpi_a2mk(x,p)
            f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] * x[i,2] + p[9] * x[i,3]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2mpi_a3(x,p)
            f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] * x[i,2] + p[9] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_Tay_a2_a2mpi_a2mk_a3(x,p)
            f = [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] * x[i,2] + p[9] * x[i,3] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end
    ##

    ## funs plot
        function model_plot_ch_a2as0(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1])
                end
            end
            return f
        end
        
        function model_plot_ch_a2as0_a3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8] * x[i,1] ^ 0.5) * x[i,1])
                end
            end
            return f
        end
        
        function model_plot_ch_a2_a2as0(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8]) * x[i,1])
                end
            end
            return f
        end
        
        function model_plot_ch_a2as0_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] + p[8] * x[i,1] * x[i,2])
                end
            end
            return f
        end
        
        function model_plot_ch_a2as0_a3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end
        
        function model_plot_ch_a2_a2as0_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end
        
        function model_plot_ch_a2as0_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] + p[8] * x[i,1] * x[i,3])
                end
            end
            return f
        end
        
        function model_plot_ch_a2as0_a3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end
        
        function model_plot_ch_a2_a2as0_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_ch_a2as1(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2as1_a3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8] * x[i,1] ^ 0.5) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as1(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8]) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2as1_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] + p[8] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_ch_a2as1_a3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as1_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_ch_a2as1_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] + p[8] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_ch_a2as1_a3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as1_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_ch_a2as2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2as2_a3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8] * x[i,1] ^ 0.5) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8]) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2as2_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] + p[8] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_ch_a2as2_a3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as2_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_ch_a2as2_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] + p[8] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_ch_a2as2_a3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as2_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_ch_a2as3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2as3_a3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8] * x[i,1] ^ 0.5) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8]) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2as3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] + p[8] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_ch_a2as3_a3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_ch_a2as3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] + p[8] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_ch_a2as3_a3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_ch_a2as4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2as4_a3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8] * x[i,1] ^ 0.5) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8]) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2as4_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] + p[8] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_ch_a2as4_a3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as4_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_ch_a2as4_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] + p[8] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_ch_a2as4_a3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as4_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_ch_a2as5(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2as5_a3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8] * x[i,1] ^ 0.5) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as5(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8]) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2as5_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] + p[8] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_ch_a2as5_a3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as5_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_ch_a2as5_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] + p[8] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_ch_a2as5_a3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as5_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end
















        function model_plot_Tay_a2as0(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1])
                end
            end
            return f
        end
        
        function model_plot_Tay_a2as0_a3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8] * x[i,1] ^ 0.5) * x[i,1])
                end
            end
            return f
        end
        
        function model_plot_Tay_a2_a2as0(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8]) * x[i,1])
                end
            end
            return f
        end
        
        function model_plot_Tay_a2as0_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] + p[8] * x[i,1] * x[i,2])
                end
            end
            return f
        end
        
        function model_plot_Tay_a2as0_a3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end
        
        function model_plot_Tay_a2_a2as0_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end
        
        function model_plot_Tay_a2as0_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0) * x[i,1] + p[8] * x[i,1] * x[i,3])
                end
            end
            return f
        end
        
        function model_plot_Tay_a2as0_a3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end
        
        function model_plot_Tay_a2_a2as0_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_0 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_Tay_a2as1(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1])
                end
            end
            return f
        end

        function model_plot_Tay_a2as1_a3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8] * x[i,1] ^ 0.5) * x[i,1])
                end
            end
            return f
        end

        function model_plot_Tay_a2_a2as1(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8]) * x[i,1])
                end
            end
            return f
        end

        function model_plot_Tay_a2as1_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] + p[8] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_Tay_a2as1_a3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_Tay_a2_a2as1_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_Tay_a2as1_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1] + p[8] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_Tay_a2as1_a3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_Tay_a2_a2as1_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_Tay_a2as2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1])
                end
            end
            return f
        end

        function model_plot_Tay_a2as2_a3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8] * x[i,1] ^ 0.5) * x[i,1])
                end
            end
            return f
        end

        function model_plot_Tay_a2_a2as2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8]) * x[i,1])
                end
            end
            return f
        end

        function model_plot_Tay_a2as2_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] + p[8] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_Tay_a2as2_a3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_Tay_a2_a2as2_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_Tay_a2as2_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2) * x[i,1] + p[8] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_Tay_a2as2_a3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_Tay_a2_a2as2_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_2 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_Tay_a2as3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1])
                end
            end
            return f
        end

        function model_plot_Tay_a2as3_a3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8] * x[i,1] ^ 0.5) * x[i,1])
                end
            end
            return f
        end

        function model_plot_Tay_a2_a2as3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8]) * x[i,1])
                end
            end
            return f
        end

        function model_plot_Tay_a2as3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] + p[8] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_Tay_a2as3_a3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_Tay_a2_a2as3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_Tay_a2as3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3) * x[i,1] + p[8] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_Tay_a2as3_a3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_Tay_a2_a2as3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_3 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_Tay_a2as4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1])
                end
            end
            return f
        end

        function model_plot_Tay_a2as4_a3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8] * x[i,1] ^ 0.5) * x[i,1])
                end
            end
            return f
        end

        function model_plot_Tay_a2_a2as4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8]) * x[i,1])
                end
            end
            return f
        end

        function model_plot_Tay_a2as4_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] + p[8] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_Tay_a2as4_a3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_Tay_a2_a2as4_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_Tay_a2as4_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4) * x[i,1] + p[8] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_Tay_a2as4_a3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_Tay_a2_a2as4_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_4 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_Tay_a2as5(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1])
                end
            end
            return f
        end

        function model_plot_Tay_a2as5_a3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8] * x[i,1] ^ 0.5) * x[i,1])
                end
            end
            return f
        end

        function model_plot_Tay_a2_a2as5(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8]) * x[i,1])
                end
            end
            return f
        end

        function model_plot_Tay_a2as5_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] + p[8] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_Tay_a2as5_a3_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_Tay_a2_a2as5_a2phi2(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,2])
                end
            end
            return f
        end

        function model_plot_Tay_a2as5_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5) * x[i,1] + p[8] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_Tay_a2as5_a3_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8] * x[i,1] ^ 0.5) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end

        function model_plot_Tay_a2_a2as5_a2phi4(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3])
                else
                    push!(f, p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_5 + p[8]) * x[i,1] + p[9] * x[i,1] * x[i,3])
                end
            end
            return f
        end



















        function model_plot_ch_a2_a2as(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2_a2as_a3(x,p)
            f = Array{uwreal,1}() 
            for i in 1:length(x[:,1])
                if x[i,1] == 0
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2))
                else
                    push!(f, p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] / (-log(sqrt(x[i,1] * value(t0_old)) * Lam)) ^ Gamma_1 + p[9] * x[i,1] ^ 0.5) * x[i,1])
                end
            end
            return f
        end

        function model_plot_ch_a2_a3(x,p)
            return [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
        end

        function model_plot_ch_a2(x,p)
            return [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7]) * x[i,1] for i in 1:length(x[:,1])]
        end

        function model_plot_ch_a2_a3_a4(x,p)
            return [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] * x[i,1] ^ 0.5 + p[9] * x[i,1]) * x[i,1] for i in 1:length(x[:,1])]
        end

        function model_plot_ch_a2_a2mpi(x,p)
            return [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
        end

        function model_plot_ch_a2_a2mpi_a2mk(x,p)
            return [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] * x[i,2] + p[9] * x[i,3]) * x[i,1] for i in 1:length(x[:,1])]
        end

        function model_plot_ch_a2_a2mpi_a3(x,p)
            return [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] * x[i,2] + p[9] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
        end

        function model_plot_ch_a2_a2mpi_a2mk_a3(x,p)
            return [p[1] * (1 + 2 * p[3] / p[2] ^ 2 * x[i,3] + (3 * p[4] - p[3]) / p[2] ^ 4 * x[i,2] ^ 2 * log(x[i,2] / p[2] ^ 2) + 4 * p[4] / p[2] ^ 4 * (x[i,3] - 0.5 * x[i,2]) ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) + (p[4] - p[3]) / p[2] ^ 4 * (4/3 * x[i,3] - x[i,2]) ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2) + 4 * p[5] / p[2] ^ 4 * x[i,3] ^ 2 + p[6] / p[2] ^ 4 * (x[i,3] - 3/2 * x[i,2]) ^ 2) + (p[7] + p[8] * x[i,2] + p[9] * x[i,3] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
        end

        function model_plot_Tay_a2_a3(x,p)
            return [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
        end

        function model_plot_Tay_a2(x,p)
            return [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7]) * x[i,1] for i in 1:length(x[:,1])]
        end

        function model_plot_Tay_a2_a3_a4(x,p)
            return [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] * x[i,1] ^ 0.5 + p[9] * x[i,1]) * x[i,1] for i in 1:length(x[:,1])]
        end

        function model_plot_Tay_a2_a2mpi(x,p)
            return [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
        end

        function model_plot_Tay_a2_a2mpi_a2mk(x,p)
            return [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] * x[i,2] + p[9] * x[i,3]) * x[i,1] for i in 1:length(x[:,1])]
        end

        function model_plot_Tay_a2_a2mpi_a3(x,p)
            return [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] * x[i,2] + p[9] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
        end

        function model_plot_Tay_a2_a2mpi_a2mk_a3(x,p)
            return [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,3] + p[5] * x[i,3] ^ 2 + p[6] * x[i,2] * x[i,3] + (p[7] + p[8] * x[i,2] + p[9] * x[i,3] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
        end
    ##

    #=
    models = [model_ch_a2, model_Tay_a2, model_ch_a2_a3, model_Tay_a2_a3, model_ch_a2_a2mpi, model_Tay_a2_a2mpi, model_ch_a2_a2mpi_a2mk, model_Tay_a2_a2mpi_a2mk, model_ch_a2_a2mpi_a3, model_Tay_a2_a2mpi_a3, model_ch_a2_a2mpi_a2mk_a3, model_Tay_a2_a2mpi_a2mk_a3]
    models_plot = [model_plot_ch_a2, model_plot_Tay_a2, model_plot_ch_a2_a3, model_plot_Tay_a2_a3, model_plot_ch_a2_a2mpi, model_plot_Tay_a2_a2mpi, model_plot_ch_a2_a2mpi_a2mk, model_plot_Tay_a2_a2mpi_a2mk, model_plot_ch_a2_a2mpi_a3, model_plot_Tay_a2_a2mpi_a3, model_plot_ch_a2_a2mpi_a2mk_a3, model_plot_Tay_a2_a2mpi_a2mk_a3]
    params = [7,7,8,8,8,8,9,9,9,9,10,10]
    =#

    models = [model_ch_a2as0, model_ch_a2as0_a3, model_ch_a2_a2as0, model_ch_a2as0_a2phi2, model_ch_a2as0_a3_a2phi2, model_ch_a2_a2as0_a2phi2, model_ch_a2as0_a2phi4, model_ch_a2as0_a3_a2phi4, model_ch_a2_a2as0_a2phi4,
              model_ch_a2as1, model_ch_a2as1_a3, model_ch_a2_a2as1, model_ch_a2as1_a2phi2, model_ch_a2as1_a3_a2phi2, model_ch_a2_a2as1_a2phi2, model_ch_a2as1_a2phi4, model_ch_a2as1_a3_a2phi4, model_ch_a2_a2as1_a2phi4,
              model_ch_a2as2, model_ch_a2as2_a3, model_ch_a2_a2as2, model_ch_a2as2_a2phi2, model_ch_a2as2_a3_a2phi2, model_ch_a2_a2as2_a2phi2, model_ch_a2as2_a2phi4, model_ch_a2as2_a3_a2phi4, model_ch_a2_a2as2_a2phi4, 
              model_ch_a2as3, model_ch_a2as3_a3, model_ch_a2_a2as3, model_ch_a2as3_a2phi2, model_ch_a2as3_a3_a2phi2, model_ch_a2_a2as3_a2phi2, model_ch_a2as3_a2phi4, model_ch_a2as3_a3_a2phi4, model_ch_a2_a2as3_a2phi4, 
              model_ch_a2as4, model_ch_a2as4_a3, model_ch_a2_a2as4, model_ch_a2as4_a2phi2, model_ch_a2as4_a3_a2phi2, model_ch_a2_a2as4_a2phi2, model_ch_a2as4_a2phi4, model_ch_a2as4_a3_a2phi4, model_ch_a2_a2as4_a2phi4, 
              model_ch_a2as5, model_ch_a2as5_a3, model_ch_a2_a2as5, model_ch_a2as5_a2phi2, model_ch_a2as5_a3_a2phi2, model_ch_a2_a2as5_a2phi2, model_ch_a2as5_a2phi4, model_ch_a2as5_a3_a2phi4, model_ch_a2_a2as5_a2phi4]

    models_plot = [model_plot_ch_a2as0, model_plot_ch_a2as0_a3, model_plot_ch_a2_a2as0, model_plot_ch_a2as0_a2phi2, model_plot_ch_a2as0_a3_a2phi2, model_plot_ch_a2_a2as0_a2phi2, model_plot_ch_a2as0_a2phi4, model_plot_ch_a2as0_a3_a2phi4, model_plot_ch_a2_a2as0_a2phi4,
                   model_plot_ch_a2as1, model_plot_ch_a2as1_a3, model_plot_ch_a2_a2as1, model_plot_ch_a2as1_a2phi2, model_plot_ch_a2as1_a3_a2phi2, model_plot_ch_a2_a2as1_a2phi2, model_plot_ch_a2as1_a2phi4, model_plot_ch_a2as1_a3_a2phi4, model_plot_ch_a2_a2as1_a2phi4, 
                   model_plot_ch_a2as2, model_plot_ch_a2as2_a3, model_plot_ch_a2_a2as2, model_plot_ch_a2as2_a2phi2, model_plot_ch_a2as2_a3_a2phi2, model_plot_ch_a2_a2as2_a2phi2, model_plot_ch_a2as2_a2phi4, model_plot_ch_a2as2_a3_a2phi4, model_plot_ch_a2_a2as2_a2phi4, 
                   model_plot_ch_a2as3, model_plot_ch_a2as3_a3, model_plot_ch_a2_a2as3, model_plot_ch_a2as3_a2phi2, model_plot_ch_a2as3_a3_a2phi2, model_plot_ch_a2_a2as3_a2phi2, model_plot_ch_a2as3_a2phi4, model_plot_ch_a2as3_a3_a2phi4, model_plot_ch_a2_a2as3_a2phi4, 
                   model_plot_ch_a2as4, model_plot_ch_a2as4_a3, model_plot_ch_a2_a2as4, model_plot_ch_a2as4_a2phi2, model_plot_ch_a2as4_a3_a2phi2, model_plot_ch_a2_a2as4_a2phi2, model_plot_ch_a2as4_a2phi4, model_plot_ch_a2as4_a3_a2phi4, model_plot_ch_a2_a2as4_a2phi4, 
                   model_plot_ch_a2as5, model_plot_ch_a2as5_a3, model_plot_ch_a2_a2as5, model_plot_ch_a2as5_a2phi2, model_plot_ch_a2as5_a3_a2phi2, model_plot_ch_a2_a2as5_a2phi2, model_plot_ch_a2as5_a2phi4, model_plot_ch_a2as5_a3_a2phi4, model_plot_ch_a2_a2as5_a2phi4]
    
    params = [7,8,8,8,9,9,8,9,9,
              7,8,8,8,9,9,8,9,9, 
              7,8,8,8,9,9,8,9,9,
              7,8,8,8,9,9,8,9,9,
              7,8,8,8,9,9,8,9,9,
              7,8,8,8,9,9,8,9,9]
    

    list_all = [
            [ens_340; ens_346; ens_355; ens_370; ens_385],
            [ens_346; ens_355; ens_370; ens_385],
            [2,3,4,5,6,8,9,10,11,12,14,15,16,17,18,19,21,22,23,25],     #nosym
            [8,9,10,11,12,14,15,16,17,18,19,21,22,23,25],           #nosym, nob340
            [1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,18,20,22,23,24,25] #mL>4.1
    ]

    plt.ioff()

    for j in 1:length(list_all)
        list = list_all[j]
        
        xx = value.([x[list,:]; x[list,:]; x[list,:]; x[list,:]])
        yy = [y[list]; x[list,1]; x[list,2]; x[list,3]]
        #xx = value.(x[list,:])
        #yy = y[list]

        #i = 1
        for i in 1:length(models)
            uprm, chi2, chi_exp, pval_aux = fit_alg(models[i], xx, yy, params[i]+3*div(length(yy),4), wpm=wpm)

            if i == 1 && j == 1 || i == 1 && j == 2 || i == 2 && j == 1 || i == 2 && j == 2
                ## plot
                    fig = figure(string(i, ",", j), (15,10))
                    subplot(131)
                    Fph = models_plot[i]([1 ./ t0 [phi2_ph for i in 1:length(t0)] phi4], uprm)
                    Fphi2 = models_plot[i]([1 ./ t0 phi2 phi4], uprm)
                    y_proj = y .- (Fphi2 .- Fph)
                    uwerr.(y_proj)
                    x_plot = [i for i in 1.08:0.01:1.35]
                    #errorbar(value.(x[ens_334,3]), value.(y_proj[ens_334]), err.(y_proj[ens_334]), err.(x[ens_334,3]), fmt="s", label=L"$\beta=3.34$", color="brown")
                    errorbar(value.(x[ens_340,3]), value.(y_proj[ens_340]), err.(y_proj[ens_340]), err.(x[ens_340,3]), fmt="s", label=L"$\beta=3.40$", color="purple")
                    aux = models_plot[i]([[1 / t0[ens_340[1]] for i in 1:length(x_plot)] [phi2_ph for i in 1:length(x_plot)] x_plot], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="purple", alpha=0.2)
                    errorbar(value.(x[ens_346,3]), value.(y_proj[ens_346]), err.(y_proj[ens_346]), err.(x[ens_346,3]), fmt="<", label=L"$\beta=3.46$", color="green")
                    aux = models_plot[i]([[1 / t0[ens_346[1]] for i in 1:length(x_plot)] [phi2_ph for i in 1:length(x_plot)] x_plot], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="green", alpha=0.2)
                    errorbar(value.(x[ens_355,3]), value.(y_proj[ens_355]), err.(y_proj[ens_355]), err.(x[ens_355,3]), fmt=">", label=L"$\beta=3.55$", color="blue")
                    aux = models_plot[i]([[1 / t0[ens_355[1]] for i in 1:length(x_plot)] [phi2_ph for i in 1:length(x_plot)] x_plot], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="blue", alpha=0.2)
                    errorbar(value.(x[ens_370,3]), value.(y_proj[ens_370]), err.(y_proj[ens_370]), err.(x[ens_370,3]), fmt="^", label=L"$\beta=3.70$", color="orange")
                    aux = models_plot[i]([[1 / t0[ens_370[1]] for i in 1:length(x_plot)] [phi2_ph for i in 1:length(x_plot)] x_plot], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="orange", alpha=0.2)
                    errorbar(value.(x[ens_385,3]), value.(y_proj[ens_385]), err.(y_proj[ens_385]), err.(x[ens_385,3]), fmt="o", label=L"$\beta=3.85$", color="red")
                    aux = models_plot[i]([[1 / t0[ens_385[1]] for i in 1:length(x_plot)] [phi2_ph for i in 1:length(x_plot)] x_plot], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="red", alpha=0.2)
                    #errorbar(value.(x[:,3]), value.(y), err.(y), err.(x[:,3]), fmt="+", label="raw")
                    legend()
                    xlabel(L"$\phi_4$")
                    ylabel(L"$\sqrt{t_0}/w_0$")
                    subplot(132)
                    Fph = models_plot[i]([1 ./ t0 phi2 [phi4_ph for i in 1:length(t0)]], uprm)
                    Fphi2 = models_plot[i]([1 ./ t0 phi2 phi4], uprm)
                    y_proj = y .- (Fphi2 .- Fph)
                    uwerr.(y_proj)
                    x_plot = [i for i in 0.07:0.01:0.8]
                    #errorbar(value.(x[ens_334,2]), value.(y_proj[ens_334]), err.(y_proj[ens_334]), err.(x[ens_334,2]), fmt="s", label=L"$\beta=3.34$", color="brown")
                    errorbar(value.(x[ens_340,2]), value.(y_proj[ens_340]), err.(y_proj[ens_340]), err.(x[ens_340,2]), fmt="s", label=L"$\beta=3.40$", color="purple")
                    aux = models_plot[i]([[1 / t0[ens_340[1]] for i in 1:length(x_plot)] x_plot [phi4_ph for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="purple", alpha=0.2)
                    errorbar(value.(x[ens_346,2]), value.(y_proj[ens_346]), err.(y_proj[ens_346]), err.(x[ens_346,2]), fmt="<", label=L"$\beta=3.46$", color="green")
                    aux = models_plot[i]([[1 / t0[ens_346[1]] for i in 1:length(x_plot)] x_plot [phi4_ph for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="green", alpha=0.2)
                    errorbar(value.(x[ens_355,2]), value.(y_proj[ens_355]), err.(y_proj[ens_355]), err.(x[ens_355,2]), fmt=">", label=L"$\beta=3.55$", color="blue")
                    aux = models_plot[i]([[1 / t0[ens_355[1]] for i in 1:length(x_plot)] x_plot [phi4_ph for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="blue", alpha=0.2)
                    errorbar(value.(x[ens_370,2]), value.(y_proj[ens_370]), err.(y_proj[ens_370]), err.(x[ens_370,2]), fmt="^", label=L"$\beta=3.70$", color="orange")
                    aux = models_plot[i]([[1 / t0[ens_370[1]] for i in 1:length(x_plot)] x_plot [phi4_ph for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="orange", alpha=0.2)
                    errorbar(value.(x[ens_385,2]), value.(y_proj[ens_385]), err.(y_proj[ens_385]), err.(x[ens_385,2]), fmt="o", label=L"$\beta=3.85$", color="red")
                    aux = models_plot[i]([[1 / t0[ens_385[1]] for i in 1:length(x_plot)] x_plot [phi4_ph for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="red", alpha=0.2)
                    #errorbar(value.(x[:,2]), value.(y), err.(y), err.(x[:,2]), fmt="+", label="raw")
                    xlabel(L"$\phi_2$")
                    subplot(133)
                    Fph = models_plot[i]([1 ./ t0 [phi2_ph for i in 1:length(t0)] [phi4_ph for i in 1:length(t0)]], uprm)
                    Fphi2 = models_plot[i]([1 ./ t0 phi2 phi4], uprm)
                    y_proj = y .- (Fphi2 .- Fph)
                    uwerr.(y_proj)
                    #errorbar(value.(x[ens_334,1]), value.(y_proj[ens_334]), err.(y_proj[ens_334]), err.(x[ens_334,1]), fmt="s", label=L"$\beta=3.34$", color="brown")
                    errorbar(value.(x[ens_340,1]), value.(y_proj[ens_340]), err.(y_proj[ens_340]), err.(x[ens_340,1]), fmt="s", label=L"$\beta=3.40$", color="purple")
                    errorbar(value.(x[ens_346,1]), value.(y_proj[ens_346]), err.(y_proj[ens_346]), err.(x[ens_346,1]), fmt="<", label=L"$\beta=3.46$", color="green")
                    errorbar(value.(x[ens_355,1]), value.(y_proj[ens_355]), err.(y_proj[ens_355]), err.(x[ens_355,1]), fmt=">", label=L"$\beta=3.55$", color="blue")
                    errorbar(value.(x[ens_370,1]), value.(y_proj[ens_370]), err.(y_proj[ens_370]), err.(x[ens_370,1]), fmt="^", label=L"$\beta=3.70$", color="orange")
                    errorbar(value.(x[ens_385,1]), value.(y_proj[ens_385]), err.(y_proj[ens_385]), err.(x[ens_385,1]), fmt="o", label=L"$\beta=3.85$", color="red")
                    #errorbar(value.(x[:,1]), value.(y), err.(y), err.(x[:,1]), fmt="+", label="raw")
                    x_plot = [i for i in 0.0:0.01:0.35]
                    aux = models_plot[i]([x_plot [phi2_ph for i in 1:length(x_plot)] [phi4_ph for i in 1:length(x_plot)]], uprm) ; uwerr.(aux)
                    v = value.(aux)
                    e = err.(aux)
                    fill_between(x_plot, v-e, v+e, color="gray", alpha=0.5)
                    xlabel(L"$a^2/t_0$")
                    tight_layout()
                    savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/ratio_i_", i, "_j_", j, ".pdf"))
                ##
            end

            #uprm, chi2, chi_exp, pval_aux = fit_alg(models[i], xx, yy, params[i], wpm=wpm)
            doff = length(yy) - (params[i]+3*div(length(yy),4))
            #doff = length(yy) - params[i]
            x_ph = [[0.0 phi2_ph phi4_ph]; [0.0 phi2_ph phi4_ph]; [0.0 phi2_ph phi4_ph]; [0.0 phi2_ph phi4_ph]]
            #x_ph = [0.0 phi2_ph phi4_ph]
            r_ph_i = models[i](x_ph, uprm)[1]; w0_ph_i = sqrt(t0_ph_ini) / r_ph_i; uwerr(r_ph_i); r_ph_i
            push!(w0_ph, w0_ph_i)
            push!(r_ph, r_ph_i)
            push!(pval_r, pval_aux)
            push!(dof_r, doff)
            push!(TIC_r, chi2 - 2 * chi_exp)
            if i == 3 && j == 2 global uprm_plot = uprm end
        end
    end
    W = exp.(-0.5 * TIC_r) / sum(exp.(-0.5 * TIC_r))
    r_ph_av = sum(r_ph .* W)
    r2_ph_av = sum(r_ph .^ 2 .* W)
    syst_r = sqrt(r2_ph_av - r_ph_av ^ 2)
    uwerr(r_ph_av)
    r_ph_av
    syst_r

    r_aux = r_ph_av + uwreal([0.0,value(syst_r)], "syst_r")
    w0_ph_av = sqrt_t0_ph_ini / r_aux; uwerr(w0_ph_av)
    sqrt_t0_av = w0_ph_ini * r_aux; uwerr(sqrt_t0_av)

    #=
    fb = BDIO_open("/home/asaez/cls_ens/results/t0_ratio.bdio", "w")
    write_uwreal(sqrt_t0_av, fb, 1)
    BDIO_close!(fb)
    =#

    uprm = uprm_plot

    uwerr.(r_ph)

    mods = [
        "[ChPT][a2]",
        "[Tay][a2]",
        "[ChPT][a2,a3]",
        "[Tay][a2,a3]",
        "[ChPT][a2,a2phi2]",
        "[Tay][a2,a2phi2]",   
        "[ChPT][a2,a2phi2,a2phi4]",
        "[Tay][a2,a2phi2,a2phi4]",
        "[ChPT][a2,a2phi2,a3]",
        "[Tay][a2,a2phi2,a3]",
        "[ChPT][a2,a2phi2,a2phi4,a3]",
        "[Tay][a2,a2phi2,a2phi4,a3]", 
        
        "[ChPT][a2][b3.40]",
        "[Tay][a2][b3.40]",
        "[ChPT][a2,a3][b3.40]",
        "[Tay][a2,a3][b3.40]",
        "[ChPT][a2,a2phi2][b3.40]",
        "[Tay][a2,a2phi2][b3.40]",   
        "[ChPT][a2,a2phi2,a2phi4][b3.40]",
        "[Tay][a2,a2phi2,a2phi4][b3.40]",
        "[ChPT][a2,a2phi2,a3][b3.40]",
        "[Tay][a2,a2phi2,a3][b3.40]",
        "[ChPT][a2,a2phi2,a2phi4,a3][b3.40]",
        "[Tay][a2,a2phi2,a2phi4,a3][b3.40]",
        
        "[ChPT][a2][sym]",
        "[Tay][a2][sym]",
        "[ChPT][a2,a3][sym]",
        "[Tay][a2,a3][sym]",
        "[ChPT][a2,a2phi2][sym]",
        "[Tay][a2,a2phi2][sym]",   
        "[ChPT][a2,a2phi2,a2phi4][sym]",
        "[Tay][a2,a2phi2,a2phi4][sym]",
        "[ChPT][a2,a2phi2,a3][sym]",
        "[Tay][a2,a2phi2,a3][sym]",
        "[ChPT][a2,a2phi2,a2phi4,a3][sym]",
        "[Tay][a2,a2phi2,a2phi4,a3][sym]",

        "[ChPT][a2][b3.40][sym]",
        "[Tay][a2][b3.40][sym]",
        "[ChPT][a2,a3][b3.40][sym]",
        "[Tay][a2,a3][b3.40][sym]",
        "[ChPT][a2,a2phi2][b3.40][sym]",
        "[Tay][a2,a2phi2][b3.40][sym]",   
        "[ChPT][a2,a2phi2,a2phi4][b3.40][sym]",
        "[Tay][a2,a2phi2,a2phi4][b3.40][sym]",
        "[ChPT][a2,a2phi2,a3][b3.40][sym]",
        "[Tay][a2,a2phi2,a3][b3.40][sym]",
        "[ChPT][a2,a2phi2,a2phi4,a3][b3.40][sym]",
        "[Tay][a2,a2phi2,a2phi4,a3][b3.40][sym]",

        "[ChPT][a2][4.1]",
        "[Tay][a2][4.1]",
        "[ChPT][a2,a3][4.1]",
        "[Tay][a2,a3][4.1]",
        "[ChPT][a2,a2phi2][4.1]",
        "[Tay][a2,a2phi2][4.1]",   
        "[ChPT][a2,a2phi2,a2phi4][4.1]",
        "[Tay][a2,a2phi2,a2phi4][4.1]",
        "[ChPT][a2,a2phi2,a3][4.1]",
        "[Tay][a2,a2phi2,a3][4.1]",
        "[ChPT][a2,a2phi2,a2phi4,a3][4.1]",
        "[Tay][a2,a2phi2,a2phi4,a3][4.1]"
    ]

    mods_as = [
        "[ChPT][a2as1]",
        "[ChPT][a2as1,a3]",
        "[ChPT][a2as2]",
        "[ChPT][a2as2,a3]",
        "[ChPT][a2as3]",
        "[ChPT][a2as3,a3]",
        "[ChPT][a2as4]",
        "[ChPT][a2as4,a3]",
        "[ChPT][a2as5]",
        "[ChPT][a2as5,a3]",
        
        "[ChPT][a2as1][b3.40]",
        "[ChPT][a2as1,a3][b3.40]",
        "[ChPT][a2as2][b3.40]",
        "[ChPT][a2as2,a3][b3.40]",
        "[ChPT][a2as3][b3.40]",
        "[ChPT][a2as3,a3][b3.40]",
        "[ChPT][a2as4][b3.40]",
        "[ChPT][a2as4,a3][b3.40]",
        "[ChPT][a2as5][b3.40]",
        "[ChPT][a2as5,a3][b3.40]",

        "[ChPT][a2as1][sym]",
        "[ChPT][a2as1,a3][sym]",
        "[ChPT][a2as2][sym]",
        "[ChPT][a2as2,a3][sym]",
        "[ChPT][a2as3][sym]",
        "[ChPT][a2as3,a3][sym]",
        "[ChPT][a2as4][sym]",
        "[ChPT][a2as4,a3][sym]",
        "[ChPT][a2as5][sym]",
        "[ChPT][a2as5,a3][sym]",

        "[ChPT][a2as1][b3.40][sym]",
        "[ChPT][a2as1,a3][b3.40][sym]",
        "[ChPT][a2as2][b3.40][sym]",
        "[ChPT][a2as2,a3][b3.40][sym]",
        "[ChPT][a2as3][b3.40][sym]",
        "[ChPT][a2as3,a3][b3.40][sym]",
        "[ChPT][a2as4][b3.40][sym]",
        "[ChPT][a2as4,a3][b3.40][sym]",
        "[ChPT][a2as5][b3.40][sym]",
        "[ChPT][a2as5,a3][b3.40][sym]",

        "[ChPT][a2as1][4.1]",
        "[ChPT][a2as1,a3][4.1]",
        "[ChPT][a2as2][4.1]",
        "[ChPT][a2as2,a3][4.1]",
        "[ChPT][a2as3][4.1]",
        "[ChPT][a2as3,a3][4.1]",
        "[ChPT][a2as4][4.1]",
        "[ChPT][a2as4,a3][4.1]",
        "[ChPT][a2as5][4.1]",
        "[ChPT][a2as5,a3][4.1]",
    ]

    fig = figure("model av")
    subplot(211)
    ax = gca() # Get the handle of the current axis
    ylabel(L"$R$")
    xx = collect(1:length(W))
    errorbar(xx, value.(r_ph), err.(r_ph), fmt="x")
    v = [value(r_ph_av) for i in 1:length(r_ph)]
    e = [sqrt(err(r_ph_av) ^ 2 + value(syst_r) ^ 2) for i in 1:length(r_ph)]
    fill_between(xx, v-e, v+e, color="deepskyblue", alpha=0.3)
    #ax[:set_ylim]([v[1]-e[1], v[1]+e[1]])
    ax[:set_ylim]([0.81, 0.83])
    plt.xticks(xx, collect(1:length(xx)))
    setp(ax.get_xticklabels(),visible=false)
    ax[:set_xlim]([0, length(xx)+1])
    subplot(212)
    ax = gca() # Get the handle of the current axis
    bar(xx, pval_r, label="pvalue")
    #plot(xx, dof_r, label="dof")
    bar(xx, W, label="weight")
    for i in div(length(xx), length(list_all)):div(length(xx), length(list_all)):length(xx)
        vlines(i,0,1,ls="--",color="red")
    end
    hlines(0.1,0,length(xx)+1,ls="--",color="red")
    ax[:set_xlim]([0, length(xx)+1])
    plt.xticks(xx, collect(1:length(xx)))
    #plt.xticks(xx, mods)
    #plt.xticks(xx, mods_as)
    xticks(rotation=90)
    #ax[:set_ylim]([0.0, maximum(W)])
    legend()
    tight_layout()

    function plot_compar_a2_r()
        errorbar(value.(x[:,1]), value.(y), err.(y), err.(x[:,1]), fmt="x", label="")
        errorbar(value.(x_Mz), value.(y_Mz), err.(y_Mz), err.(x_Mz), fmt="+", label="Mainz")
        xlabel(L"$a^2/t_0$")
        ylabel(L"$\sqrt{t_0}/w_0$")
        legend()
    end

    function plot_compar_phi2_r()
        errorbar(value.(x[:,2]), value.(y), err.(y), err.(x[:,2]), fmt="x", label="")
        errorbar(value.(phi2_Mz), value.(y_Mz), err.(y_Mz), err.(phi2_Mz), fmt="+", label="Mainz")
        xlabel(L"$\phi_2$")
        ylabel(L"$\sqrt{t_0}/w_0$")
        legend()
    end

    function plot_compar_phi4_r()
        errorbar(value.(x[:,3]), value.(y), err.(y), err.(x[:,3]), fmt="x", label="")
        errorbar(value.(phi4_Mz), value.(y_Mz), err.(y_Mz), err.(phi4_Mz), fmt="+", label="Mainz")
        xlabel(L"$\phi_4$")
        ylabel(L"$\sqrt{t_0}/w_0$")
        legend()
    end

    function proj_phi2_r()
        Fph = models_plot[i]([1 ./ t0 [phi2_ph for i in 1:length(t0)] phi4], uprm)
        Fphi2 = models_plot[i]([1 ./ t0 phi2 phi4], uprm)
        y_proj = y .- (Fphi2 .- Fph)
        uwerr.(y_proj)

        x_plot = [i for i in 1.08:0.01:1.35]
        #errorbar(value.(x[ens_334,3]), value.(y_proj[ens_334]), err.(y_proj[ens_334]), err.(x[ens_334,3]), fmt="s", label=L"$\beta=3.34$", color="brown")
        errorbar(value.(x[ens_340,3]), value.(y_proj[ens_340]), err.(y_proj[ens_340]), err.(x[ens_340,3]), fmt="s", label=L"$\beta=3.40$", color="purple")
        aux = models_plot[i]([[1 / t0[ens_340[1]] for i in 1:length(x_plot)] [phi2_ph for i in 1:length(x_plot)] x_plot], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="purple", alpha=0.2)
        errorbar(value.(x[ens_346,3]), value.(y_proj[ens_346]), err.(y_proj[ens_346]), err.(x[ens_346,3]), fmt="<", label=L"$\beta=3.46$", color="green")
        aux = models_plot[i]([[1 / t0[ens_346[1]] for i in 1:length(x_plot)] [phi2_ph for i in 1:length(x_plot)] x_plot], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="green", alpha=0.2)
        errorbar(value.(x[ens_355,3]), value.(y_proj[ens_355]), err.(y_proj[ens_355]), err.(x[ens_355,3]), fmt=">", label=L"$\beta=3.55$", color="blue")
        aux = models_plot[i]([[1 / t0[ens_355[1]] for i in 1:length(x_plot)] [phi2_ph for i in 1:length(x_plot)] x_plot], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="blue", alpha=0.2)
        errorbar(value.(x[ens_370,3]), value.(y_proj[ens_370]), err.(y_proj[ens_370]), err.(x[ens_370,3]), fmt="^", label=L"$\beta=3.70$", color="orange")
        aux = models_plot[i]([[1 / t0[ens_370[1]] for i in 1:length(x_plot)] [phi2_ph for i in 1:length(x_plot)] x_plot], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="orange", alpha=0.2)
        errorbar(value.(x[ens_385,3]), value.(y_proj[ens_385]), err.(y_proj[ens_385]), err.(x[ens_385,3]), fmt="o", label=L"$\beta=3.85$", color="red")
        aux = models_plot[i]([[1 / t0[ens_385[1]] for i in 1:length(x_plot)] [phi2_ph for i in 1:length(x_plot)] x_plot], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="red", alpha=0.2)

        #errorbar(value.(x[:,3]), value.(y), err.(y), err.(x[:,3]), fmt="+", label="raw")

        legend()
        xlabel(L"$\phi_4$")
        ylabel(L"$\sqrt{t_0}/w_0$")
        tight_layout()
    end

    function proj_phi4_r()
        Fph = models_plot[i]([1 ./ t0 phi2 [phi4_ph for i in 1:length(t0)]], uprm)
        Fphi2 = models_plot[i]([1 ./ t0 phi2 phi4], uprm)
        y_proj = y .- (Fphi2 .- Fph)
        uwerr.(y_proj)

        x_plot = [i for i in 0.07:0.01:0.8]
        #errorbar(value.(x[ens_334,2]), value.(y_proj[ens_334]), err.(y_proj[ens_334]), err.(x[ens_334,2]), fmt="s", label=L"$\beta=3.34$", color="brown")
        errorbar(value.(x[ens_340,2]), value.(y_proj[ens_340]), err.(y_proj[ens_340]), err.(x[ens_340,2]), fmt="s", label=L"$\beta=3.40$", color="purple")
        aux = models_plot[i]([[1 / t0[ens_340[1]] for i in 1:length(x_plot)] x_plot [phi4_ph for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="purple", alpha=0.2)
        errorbar(value.(x[ens_346,2]), value.(y_proj[ens_346]), err.(y_proj[ens_346]), err.(x[ens_346,2]), fmt="<", label=L"$\beta=3.46$", color="green")
        aux = models_plot[i]([[1 / t0[ens_346[1]] for i in 1:length(x_plot)] x_plot [phi4_ph for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="green", alpha=0.2)
        errorbar(value.(x[ens_355,2]), value.(y_proj[ens_355]), err.(y_proj[ens_355]), err.(x[ens_355,2]), fmt=">", label=L"$\beta=3.55$", color="blue")
        aux = models_plot[i]([[1 / t0[ens_355[1]] for i in 1:length(x_plot)] x_plot [phi4_ph for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="blue", alpha=0.2)
        errorbar(value.(x[ens_370,2]), value.(y_proj[ens_370]), err.(y_proj[ens_370]), err.(x[ens_370,2]), fmt="^", label=L"$\beta=3.70$", color="orange")
        aux = models_plot[i]([[1 / t0[ens_370[1]] for i in 1:length(x_plot)] x_plot [phi4_ph for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="orange", alpha=0.2)
        errorbar(value.(x[ens_385,2]), value.(y_proj[ens_385]), err.(y_proj[ens_385]), err.(x[ens_385,2]), fmt="o", label=L"$\beta=3.85$", color="red")
        aux = models_plot[i]([[1 / t0[ens_385[1]] for i in 1:length(x_plot)] x_plot [phi4_ph for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="red", alpha=0.2)
        #errorbar(value.(x[:,2]), value.(y), err.(y), err.(x[:,2]), fmt="+", label="raw")
        legend()
        xlabel(L"$\phi_2$")
        ylabel(L"$\sqrt{t_0}/w_0$")
        tight_layout()
    end

    function proj_phi24_r()
        Fph = models_plot[i]([1 ./ t0 [phi2_ph for i in 1:length(t0)] [phi4_ph for i in 1:length(t0)]], uprm)
        Fphi2 = models_plot[i]([1 ./ t0 phi2 phi4], uprm)
        y_proj = y .- (Fphi2 .- Fph)
        uwerr.(y_proj)

        #errorbar(value.(x[ens_334,1]), value.(y_proj[ens_334]), err.(y_proj[ens_334]), err.(x[ens_334,1]), fmt="s", label=L"$\beta=3.34$", color="brown")
        errorbar(value.(x[ens_340,1]), value.(y_proj[ens_340]), err.(y_proj[ens_340]), err.(x[ens_340,1]), fmt="s", label=L"$\beta=3.40$", color="purple")
        errorbar(value.(x[ens_346,1]), value.(y_proj[ens_346]), err.(y_proj[ens_346]), err.(x[ens_346,1]), fmt="<", label=L"$\beta=3.46$", color="green")
        errorbar(value.(x[ens_355,1]), value.(y_proj[ens_355]), err.(y_proj[ens_355]), err.(x[ens_355,1]), fmt=">", label=L"$\beta=3.55$", color="blue")
        errorbar(value.(x[ens_370,1]), value.(y_proj[ens_370]), err.(y_proj[ens_370]), err.(x[ens_370,1]), fmt="^", label=L"$\beta=3.70$", color="orange")
        errorbar(value.(x[ens_385,1]), value.(y_proj[ens_385]), err.(y_proj[ens_385]), err.(x[ens_385,1]), fmt="o", label=L"$\beta=3.85$", color="red")
        #errorbar(value.(x[:,1]), value.(y), err.(y), err.(x[:,1]), fmt="+", label="raw")
        legend()
        x_plot = [i for i in 0.0:0.01:0.35]
        aux = models_plot[i]([x_plot [phi2_ph for i in 1:length(x_plot)] [phi4_ph for i in 1:length(x_plot)]], uprm) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot, v-e, v+e, color="gray", alpha=0.5)
        legend()
        xlabel(L"$a^2/t_0$")
        ylabel(L"$\sqrt{t_0}/w_0$")
        tight_layout()
    end

#===========================================#

#=========== w0 Omega =================#

    read_data("omega")

    L = [ens_db_Mainz[ens_av_m[i]][1] for i in 1:length(ens_av_m)]
    x = [1 ./ w0 wphi2 wphi4 mpi .* L]; uwerr.(x)
    y = sqrt.(w0) .* mO; uwerr.(y)

    TIC_w0 = Array{Float64,1}()
    pval_w0 = Array{Float64,1}()
    dof_w0 = Array{Float64,1}()
    w0_ph = Array{uwreal,1}()

    ## funs
        function model_ch_a2(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a3(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a3_a4(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,1] ^ 0.5 + p[10] * x[i,1]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2mpi(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2mpi_a2mk(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2] + p[10] * x[i,3]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2mpi_a3(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2mpi_a2mk_a3(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2] + p[10] * x[i,3] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_a2(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_a2_a3(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_a2_a3_a4(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,1] ^ 0.5 + p[11] * x[i,1]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_a2_a2mpi(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_a2_a2mpi_a2mk(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_a2_a2mpi_a3(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_a2_a2mpi_a2mk_a3(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3] + p[12] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+12] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_fve_a2(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_fve_a2_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_fve_a2_a3_a4(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5 + p[11] * x[i,1]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_fve_a2_a2mpi(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_fve_a2_a2mpi_a2mk(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_fve_a2_a2mpi_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_fve_a2_a2mpi_a2mk_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3] + p[12] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+12] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end
        
        function model_ch_fve_a2(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_fve_a2_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_fve_a2_a3_a4(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5 + p[11] * x[i,1]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_fve_a2_a2mpi(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_fve_a2_a2mpi_a2mk(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_fve_a2_a2mpi_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_fve_a2_a2mpi_a2mk_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3] + p[12] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+12] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end
    ##

    ## funs plot
        function model_plot_ch_a2(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_a2_a3(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_a2_a3_a4(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,1] ^ 0.5 + p[10] * x[i,1]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_a2_a2mpi(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_a2_a2mpi_a2mk(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2] + p[10] * x[i,3]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_a2_a2mpi_a3(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_a2_a2mpi_a2mk_a3(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2] + p[10] * x[i,3] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_a2(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_a2_a3(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_a2_a3_a4(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,1] ^ 0.5 + p[11] * x[i,1]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_a2_a2mpi(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_a2_a2mpi_a2mk(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_a2_a2mpi_a3(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_a2_a2mpi_a2mk_a3(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3] + p[12] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+12] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_fve_a2(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_fve_a2_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_fve_a2_a3_a4(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5 + p[11] * x[i,1]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_fve_a2_a2mpi(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_fve_a2_a2mpi_a2mk(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_fve_a2_a2mpi_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_fve_a2_a2mpi_a2mk_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3] + p[12] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+12] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_fve_a2(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_fve_a2_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_fve_a2_a3_a4(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5 + p[11] * x[i,1]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_fve_a2_a2mpi(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_fve_a2_a2mpi_a2mk(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_fve_a2_a2mpi_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_fve_a2_a2mpi_a2mk_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3] + p[12] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+12] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end
    ##

    models = [model_ch_a2, model_ch_a2_a3, model_ch_a2_a2mpi, model_ch_a2_a2mpi_a2mk, model_ch_a2_a2mpi_a3, model_ch_a2_a2mpi_a2mk_a3,
              model_ch_log_a2, model_ch_log_a2_a3, model_ch_log_a2_a2mpi, model_ch_log_a2_a2mpi_a2mk, model_ch_log_a2_a2mpi_a3, model_ch_log_a2_a2mpi_a2mk_a3,
              model_ch_log_fve_a2, model_ch_log_fve_a2_a3, model_ch_log_fve_a2_a2mpi, model_ch_log_fve_a2_a2mpi_a2mk, model_ch_log_fve_a2_a2mpi_a3, model_ch_log_fve_a2_a2mpi_a2mk_a3,
              model_ch_fve_a2, model_ch_fve_a2_a3, model_ch_fve_a2_a2mpi, model_ch_fve_a2_a2mpi_a2mk, model_ch_fve_a2_a2mpi_a3, model_ch_fve_a2_a2mpi_a2mk_a3
    ]
    models_plot = [model_plot_ch_a2, model_plot_ch_a2_a3, model_plot_ch_a2_a2mpi, model_plot_ch_a2_a2mpi_a2mk, model_plot_ch_a2_a2mpi_a3, model_plot_ch_a2_a2mpi_a2mk_a3,
                   model_plot_ch_log_a2, model_plot_ch_log_a2_a3, model_plot_ch_log_a2_a2mpi, model_plot_ch_log_a2_a2mpi_a2mk, model_plot_ch_log_a2_a2mpi_a3, model_plot_ch_log_a2_a2mpi_a2mk_a3,
                   model_plot_ch_log_fve_a2, model_plot_ch_log_fve_a2_a3, model_plot_ch_log_fve_a2_a2mpi, model_plot_ch_log_fve_a2_a2mpi_a2mk, model_plot_ch_log_fve_a2_a2mpi_a3, model_plot_ch_log_fve_a2_a2mpi_a2mk_a3,
                   model_plot_ch_fve_a2, model_plot_ch_fve_a2_a3, model_plot_ch_fve_a2_a2mpi, model_plot_ch_fve_a2_a2mpi_a2mk, model_plot_ch_fve_a2_a2mpi_a3, model_plot_ch_fve_a2_a2mpi_a2mk_a3
    ]
    params = [7,9,9,10,10,11,9,10,10,11,11,12,9,10,10,11,11,12,9,10,10,11,11,12]

    list_all = [
            [ens_340; ens_346; ens_355; ens_370; ens_385],
            [ens_346; ens_355; ens_370; ens_385],
            [2,3,4,6,7,8,9,11,12,13,14,16,17,19],     #nosym
            [1,2,3,4,5,6,7,8,10,11,12,13,15,17,18,19] #mL>4.1
    ]

    for j in 1:length(list_all)
        list = list_all[j]
        xx = value.([x[list,:]; x[list,:]; x[list,:]; x[list,:]])
        yy = [y[list]; x[list,1]; x[list,2]; x[list,3]]

        #i = 1
        for i in 1:length(models)
            uprm, chi2, chi_exp, pval_aux = fit_alg(models[i], xx, yy, params[i]+3*div(length(yy),4), wpm=wpm)
            doff = length(yy) - (params[i]+3*div(length(yy),4))

            function fun(x,p)
                a = [[0.0,p[1],p[2],1000]'; [0.0,p[1],p[2],1000]'; [0.0,p[1],p[2],1000]'; [0.0,p[1],p[2],1000]']
                f = models[i](a,x)[1] .^ 2 ./ p'[:,1]
                g = models[i](a,x)[1] .^ 2 ./ p'[:,2]
                h = [p[2+k] for k in 1:div(length(x),2)]
                return [f;g;h]
            end

            #a = fit_alg(fun, [value.(uprm);value.(uprm)], [ratio1_ph;ratio2_ph;uprm], 2+length(uprm), [value(wphi2_ph),value(wphi4_ph)])
            #x_ph = [[0.0 a[1][1] a[1][2] 1000]; [0.0 a[1][1] a[1][2] 1000]; [0.0 a[1][1] a[1][2] 1000]; [0.0 a[1][1] a[1][2] 1000]]
            x_ph = [[0.0 wphi2_ph wphi4_ph 1000]; [0.0 wphi2_ph wphi4_ph 1000]; [0.0 wphi2_ph wphi4_ph 1000]; [0.0 wphi2_ph wphi4_ph 1000]]
            w0mO_ph = models[i](x_ph, uprm)[1]; w0_ph_i = w0mO_ph / mO_ph * hc; uwerr(w0_ph_i); w0_ph_i
                        
            push!(w0_ph, w0_ph_i)
            push!(pval_w0, pval_aux)
            push!(dof_w0, doff)
            push!(TIC_w0, chi2 - 2 * chi_exp)

            if i == 3 && j == 1 global uprm_plot = uprm end
        end
    end
    W = exp.(-0.5 * TIC_w0) / sum(exp.(-0.5 * TIC_w0))
    w0_ph_av = sum(w0_ph .* W)
    w02_ph_av = sum(w0_ph .^ 2 .* W)
    syst_w0 = sqrt(w02_ph_av - w0_ph_av ^ 2)
    uwerr(w0_ph_av)
    w0_ph_av
    syst_w0
    #w0_ph_ini = deepcopy(w0_ph)

    uprm = uprm_plot

    function proj_psi2_w0mO()
        Fph = models_plot[i]([1 ./ w0 [wphi2_ph for i in 1:length(w0)] wphi4 [1000 for i in 1:length(t0)]], uprm)
        Fphi2 = models_plot[i]([1 ./ w0 wphi2 wphi4 value.(x[:,4])], uprm)
        y_proj = y .- (Fphi2 .- Fph)
        uwerr.(y_proj)

        x_plot = [i for i in 1.43:0.01:1.6]
        #errorbar(value.(x[ens_334,3]), value.(y_proj[ens_334]), err.(y_proj[ens_334]), err.(x[ens_334,3]), fmt="s", label=L"$\beta=3.34$", color="brown")
        errorbar(value.(x[ens_340,3]), value.(y_proj[ens_340]), err.(y_proj[ens_340]), err.(x[ens_340,3]), fmt="s", label=L"$\beta=3.40$", color="purple")
        aux = models_plot[i]([[1 / w0[ens_340[1]] for i in 1:length(x_plot)] [wphi2_ph for i in 1:length(x_plot)] x_plot [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="purple", alpha=0.2)
        errorbar(value.(x[ens_346,3]), value.(y_proj[ens_346]), err.(y_proj[ens_346]), err.(x[ens_346,3]), fmt="<", label=L"$\beta=3.46$", color="green")
        aux = models_plot[i]([[1 / w0[ens_346[1]] for i in 1:length(x_plot)] [wphi2_ph for i in 1:length(x_plot)] x_plot [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="green", alpha=0.2)
        errorbar(value.(x[ens_355,3]), value.(y_proj[ens_355]), err.(y_proj[ens_355]), err.(x[ens_355,3]), fmt=">", label=L"$\beta=3.55$", color="blue")
        aux = models_plot[i]([[1 / w0[ens_355[1]] for i in 1:length(x_plot)] [wphi2_ph for i in 1:length(x_plot)] x_plot [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="blue", alpha=0.2)
        errorbar(value.(x[ens_370,3]), value.(y_proj[ens_370]), err.(y_proj[ens_370]), err.(x[ens_370,3]), fmt="^", label=L"$\beta=3.70$", color="orange")
        aux = models_plot[i]([[1 / w0[ens_370[1]] for i in 1:length(x_plot)] [wphi2_ph for i in 1:length(x_plot)] x_plot [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="orange", alpha=0.2)
        errorbar(value.(x[ens_385,3]), value.(y_proj[ens_385]), err.(y_proj[ens_385]), err.(x[ens_385,3]), fmt="o", label=L"$\beta=3.85$", color="red")
        aux = models_plot[i]([[1 / w0[ens_385[1]] for i in 1:length(x_plot)] [wphi2_ph for i in 1:length(x_plot)] x_plot [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="red", alpha=0.2)

        #errorbar(value.(x[:,3]), value.(y), err.(y), err.(x[:,3]), fmt="+", label="raw")

        legend()
        xlabel(L"$\psi_4$")
        ylabel(L"$w_0m_{\Omega}$")
        tight_layout()
    end

    function proj_psi4_w0mO()
        Fph = models_plot[i]([1 ./ w0 wphi2 [wphi4_ph for i in 1:length(w0)] [1000 for i in 1:length(t0)]], uprm)
        Fphi2 = models_plot[i]([1 ./ w0 wphi2 wphi4 value.(x[:,4])], uprm)
        y_proj = y .- (Fphi2 .- Fph)
        uwerr.(y_proj)

        x_plot = [i for i in 0.1:0.01:1.1]
        #errorbar(value.(x[ens_334,2]), value.(y_proj[ens_334]), err.(y_proj[ens_334]), err.(x[ens_334,2]), fmt="s", label=L"$\beta=3.34$", color="brown")
        errorbar(value.(x[ens_340,2]), value.(y_proj[ens_340]), err.(y_proj[ens_340]), err.(x[ens_340,2]), fmt="s", label=L"$\beta=3.40$", color="purple")
        aux = models_plot[i]([[1 / w0[ens_340[1]] for i in 1:length(x_plot)] x_plot [wphi4_ph for i in 1:length(x_plot)] [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="purple", alpha=0.2)
        errorbar(value.(x[ens_346,2]), value.(y_proj[ens_346]), err.(y_proj[ens_346]), err.(x[ens_346,2]), fmt="<", label=L"$\beta=3.46$", color="green")
        aux = models_plot[i]([[1 / w0[ens_346[1]] for i in 1:length(x_plot)] x_plot [wphi4_ph for i in 1:length(x_plot)] [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="green", alpha=0.2)
        errorbar(value.(x[ens_355,2]), value.(y_proj[ens_355]), err.(y_proj[ens_355]), err.(x[ens_355,2]), fmt=">", label=L"$\beta=3.55$", color="blue")
        aux = models_plot[i]([[1 / w0[ens_355[1]] for i in 1:length(x_plot)] x_plot [wphi4_ph for i in 1:length(x_plot)] [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="blue", alpha=0.2)
        errorbar(value.(x[ens_370,2]), value.(y_proj[ens_370]), err.(y_proj[ens_370]), err.(x[ens_370,2]), fmt="^", label=L"$\beta=3.70$", color="orange")
        aux = models_plot[i]([[1 / w0[ens_370[1]] for i in 1:length(x_plot)] x_plot [wphi4_ph for i in 1:length(x_plot)] [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="orange", alpha=0.2)
        errorbar(value.(x[ens_385,2]), value.(y_proj[ens_385]), err.(y_proj[ens_385]), err.(x[ens_385,2]), fmt="o", label=L"$\beta=3.85$", color="red")
        aux = models_plot[i]([[1 / w0[ens_385[1]] for i in 1:length(x_plot)] x_plot [wphi4_ph for i in 1:length(x_plot)] [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="red", alpha=0.2)
        #errorbar(value.(x[:,2]), value.(y), err.(y), err.(x[:,2]), fmt="+", label="raw")
        legend()
        xlabel(L"$\psi_2$")
        ylabel(L"$w_0m_{\Omega}$")
        tight_layout()
    end

    function proj_psi24_w0mO()
        Fph = models_plot[i]([1 ./ w0 [wphi2_ph for i in 1:length(t0)] [wphi4_ph for i in 1:length(w0)] [1000 for i in 1:length(t0)]], uprm)
        Fphi2 = models_plot[i]([1 ./ w0 wphi2 wphi4 value.(x[:,4])], uprm)
        y_proj = y .- (Fphi2 .- Fph)
        uwerr.(y_proj)

        #errorbar(value.(x[ens_334,1]), value.(y_proj[ens_334]), err.(y_proj[ens_334]), err.(x[ens_334,1]), fmt="s", label=L"$\beta=3.34$", color="brown")
        errorbar(value.(x[ens_340,1]), value.(y_proj[ens_340]), err.(y_proj[ens_340]), err.(x[ens_340,1]), fmt="s", label=L"$\beta=3.40$", color="purple")
        errorbar(value.(x[ens_346,1]), value.(y_proj[ens_346]), err.(y_proj[ens_346]), err.(x[ens_346,1]), fmt="<", label=L"$\beta=3.46$", color="green")
        errorbar(value.(x[ens_355,1]), value.(y_proj[ens_355]), err.(y_proj[ens_355]), err.(x[ens_355,1]), fmt=">", label=L"$\beta=3.55$", color="blue")
        errorbar(value.(x[ens_370,1]), value.(y_proj[ens_370]), err.(y_proj[ens_370]), err.(x[ens_370,1]), fmt="^", label=L"$\beta=3.70$", color="orange")
        errorbar(value.(x[ens_385,1]), value.(y_proj[ens_385]), err.(y_proj[ens_385]), err.(x[ens_385,1]), fmt="o", label=L"$\beta=3.85$", color="red")
        #errorbar(value.(x[:,1]), value.(y), err.(y), err.(x[:,1]), fmt="+", label="raw")
        legend()
        x_plot = [i for i in 0.0:0.01:0.3]
        aux = models_plot[i]([x_plot [wphi2_ph for i in 1:length(x_plot)] [wphi4_ph for i in 1:length(x_plot)] [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot, v-e, v+e, color="gray", alpha=0.5)
        legend()
        xlabel(L"$a^2/t_0$")
        ylabel(L"$w_0m_{\Omega}$")
        tight_layout()
    end

    mods = [
        "[no log][no fve][a2]",
        "[no log][no fve][a2, a3]",
        "[no log][no fve][a2, a2phi2]",
        "[no log][no fve][a2, a2phi2, a2phi4]",
        "[no log][no fve][a2, a2phi2, a3]",
        "[no log][no fve][a2, a2phi2, a2phi4, a3]",
        "[log][no fve][a2]",
        "[log][no fve][a2, a3]",
        "[log][no fve][a2, a2phi2]",
        "[log][no fve][a2, a2phi2, a2phi4]",
        "[log][no fve][a2, a2phi2, a3]",
        "[log][no fve][a2, a2phi2, a2phi4, a3]",
        "[log][fve][a2]",
        "[log][fve][a2, a3]",
        "[log][fve][a2, a2phi2]",
        "[log][fve][a2, a2phi2, a2phi4]",
        "[log][fve][a2, a2phi2, a3]",
        "[log][fve][a2, a2phi2, a2phi4, a3]",
        "[no log][fve][a2]",
        "[no log][fve][a2, a3]",
        "[no log][fve][a2, a2phi2]",
        "[no log][fve][a2, a2phi2, a2phi4]",
        "[no log][fve][a2, a2phi2, a3]",
        "[no log][fve][a2, a2phi2, a2phi4, a3]",

        "[no log][no fve][a2][b3.40]",
        "[no log][no fve][a2, a3][b3.40]",
        "[no log][no fve][a2, a2phi2][b3.40]",
        "[no log][no fve][a2, a2phi2, a2phi4][b3.40]",
        "[no log][no fve][a2, a2phi2, a3][b3.40]",
        "[no log][no fve][a2, a2phi2, a2phi4, a3][b3.40]",
        "[log][no fve][a2][b3.40]",
        "[log][no fve][a2, a3][b3.40]",
        "[log][no fve][a2, a2phi2][b3.40]",
        "[log][no fve][a2, a2phi2, a2phi4][b3.40]",
        "[log][no fve][a2, a2phi2, a3][b3.40]",
        "[log][no fve][a2, a2phi2, a2phi4, a3][b3.40]",
        "[log][fve][a2][b3.40]",
        "[log][fve][a2, a3][b3.40]",
        "[log][fve][a2, a2phi2][b3.40]",
        "[log][fve][a2, a2phi2, a2phi4][b3.40]",
        "[log][fve][a2, a2phi2, a3][b3.40]",
        "[log][fve][a2, a2phi2, a2phi4, a3][b3.40]",
        "[no log][fve][a2][b3.40]",
        "[no log][fve][a2, a3][b3.40]",
        "[no log][fve][a2, a2phi2][b3.40]",
        "[no log][fve][a2, a2phi2, a2phi4][b3.40]",
        "[no log][fve][a2, a2phi2, a3][b3.40]",
        "[no log][fve][a2, a2phi2, a2phi4, a3][b3.40]",

        "[no log][no fve][a2][sym]",
        "[no log][no fve][a2, a3][sym]",
        "[no log][no fve][a2, a2phi2][sym]",
        "[no log][no fve][a2, a2phi2, a2phi4][sym]",
        "[no log][no fve][a2, a2phi2, a3][sym]",
        "[no log][no fve][a2, a2phi2, a2phi4, a3][sym]",
        "[log][no fve][a2][sym]",
        "[log][no fve][a2, a3][sym]",
        "[log][no fve][a2, a2phi2][sym]",
        "[log][no fve][a2, a2phi2, a2phi4][sym]",
        "[log][no fve][a2, a2phi2, a3][sym]",
        "[log][no fve][a2, a2phi2, a2phi4, a3][sym]",
        "[log][fve][a2][sym]",
        "[log][fve][a2, a3][sym]",
        "[log][fve][a2, a2phi2][sym]",
        "[log][fve][a2, a2phi2, a2phi4][sym]",
        "[log][fve][a2, a2phi2, a3][sym]",
        "[log][fve][a2, a2phi2, a2phi4, a3][sym]",
        "[no log][fve][a2][sym]",
        "[no log][fve][a2, a3][sym]",
        "[no log][fve][a2, a2phi2][sym]",
        "[no log][fve][a2, a2phi2, a2phi4][sym]",
        "[no log][fve][a2, a2phi2, a3][sym]",
        "[no log][fve][a2, a2phi2, a2phi4, a3][sym]",

        "[no log][no fve][a2][4.1]",
        "[no log][no fve][a2, a3][4.1]",
        "[no log][no fve][a2, a2phi2][4.1]",
        "[no log][no fve][a2, a2phi2, a2phi4][4.1]",
        "[no log][no fve][a2, a2phi2, a3][4.1]",
        "[no log][no fve][a2, a2phi2, a2phi4, a3][4.1]",
        "[log][no fve][a2][4.1]",
        "[log][no fve][a2, a3][4.1]",
        "[log][no fve][a2, a2phi2][4.1]",
        "[log][no fve][a2, a2phi2, a2phi4][4.1]",
        "[log][no fve][a2, a2phi2, a3][4.1]",
        "[log][no fve][a2, a2phi2, a2phi4, a3][4.1]",
        "[log][fve][a2][4.1]",
        "[log][fve][a2, a3][4.1]",
        "[log][fve][a2, a2phi2][4.1]",
        "[log][fve][a2, a2phi2, a2phi4][4.1]",
        "[log][fve][a2, a2phi2, a3][4.1]",
        "[log][fve][a2, a2phi2, a2phi4, a3][4.1]",
        "[no log][fve][a2][4.1]",
        "[no log][fve][a2, a3][4.1]",
        "[no log][fve][a2, a2phi2][4.1]",
        "[no log][fve][a2, a2phi2, a2phi4][4.1]",
        "[no log][fve][a2, a2phi2, a3][4.1]",
        "[no log][fve][a2, a2phi2, a2phi4, a3][4.1]"
    ]

    fig = figure("model av")
    subplot(211)
    ax = gca() # Get the handle of the current axis
    ylabel(L"$w_0\; [fm]$")
    xx = collect(1:length(W))
    errorbar(xx, value.(w0_ph), err.(w0_ph), fmt="x")
    v = [value(w0_ph_av) for i in 1:length(w0_ph)]
    e = [sqrt(err(w0_ph_av) ^ 2 + value(syst_w0) ^ 2) for i in 1:length(w0_ph)]
    fill_between(xx, v-e, v+e, color="deepskyblue", alpha=0.3)
    #ax[:set_ylim]([v[1]-e[1], v[1]+e[1]])
    ax[:set_ylim]([0.16, 0.18])
    plt.xticks(xx, collect(1:length(xx)))
    setp(ax.get_xticklabels(),visible=false)
    ax[:set_xlim]([0, length(xx)+1])
    subplot(212)
    ax = gca() # Get the handle of the current axis
    bar(xx, pval_w0, label="pvalue")
    #plot(xx, dof_w0, label="dof")
    bar(xx, W, label="weight")
    for i in div(length(xx), length(list_all)):div(length(xx), length(list_all)):length(xx)
        vlines(i,0,1,ls="--",color="red")
    end
    ax[:set_xlim]([0, length(xx)+1])
    #plt.xticks(xx, collect(1:length(xx)))
    plt.xticks(xx, mods)
    xticks(rotation=90)
    #ax[:set_ylim]([0.0, maximum(W)])
    legend()
    tight_layout()

#===========================================#

#=========== t0 Omega =================#

    read_data("omega")

    L = [ens_db_Mainz[ens_av_m[i]][1] for i in 1:length(ens_av_m)]
    x = [1 ./ t0 phi2 phi4 mpi .* L]; uwerr.(x)
    y = sqrt.(t0) .* mO; uwerr.(y)

    TIC_t0 = Array{Float64,1}()
    pval_t0 = Array{Float64,1}()
    dof_t0 = Array{Float64,1}()
    t0_ph = Array{uwreal,1}()

    ## funs
        function model_ch_a2(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a3(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a3_a4(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,1] ^ 0.5 + p[10] * x[i,1]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2mpi(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2mpi_a2mk(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2] + p[10] * x[i,3]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2mpi_a3(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_a2_a2mpi_a2mk_a3(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2] + p[10] * x[i,3] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_a2(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_a2_a3(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_a2_a3_a4(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,1] ^ 0.5 + p[11] * x[i,1]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_a2_a2mpi(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_a2_a2mpi_a2mk(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_a2_a2mpi_a3(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_a2_a2mpi_a2mk_a3(x,p)
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3] + p[12] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+12] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_fve_a2(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_fve_a2_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_fve_a2_a3_a4(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5 + p[11] * x[i,1]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_fve_a2_a2mpi(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_fve_a2_a2mpi_a2mk(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_fve_a2_a2mpi_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_log_fve_a2_a2mpi_a2mk_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3] + p[12] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+12] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end
        
        function model_ch_fve_a2(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_fve_a2_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_fve_a2_a3_a4(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5 + p[11] * x[i,1]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_fve_a2_a2mpi(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_fve_a2_a2mpi_a2mk(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3]) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_fve_a2_a2mpi_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end

        function model_ch_fve_a2_a2mpi_a2mk_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3] + p[12] * x[i,1] ^ 0.5) * x[i,1] for i in 1:div(length(x[:,1]),4)]
            g = [p[i+12] for i in 1:3*div(length(x[:,1]),4)]
            return [f;g]
        end
    ##

    ## funs plot
        function model_plot_ch_a2(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+7] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_a2_a3(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_a2_a3_a4(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,1] ^ 0.5 + p[10] * x[i,1]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_a2_a2mpi(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_a2_a2mpi_a2mk(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2] + p[10] * x[i,3]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_a2_a2mpi_a3(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_a2_a2mpi_a2mk_a3(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4  + (p[8] + p[9] * x[i,2] + p[10] * x[i,3] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_a2(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_a2_a3(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_a2_a3_a4(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,1] ^ 0.5 + p[11] * x[i,1]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_a2_a2mpi(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_a2_a2mpi_a2mk(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_a2_a2mpi_a3(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_a2_a2mpi_a2mk_a3(x,p)
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * log(x[i,2] / p[1] ^ 2) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3] + p[12] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+12] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_fve_a2(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_fve_a2_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_fve_a2_a3_a4(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5 + p[11] * x[i,1]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_fve_a2_a2mpi(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_fve_a2_a2mpi_a2mk(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_fve_a2_a2mpi_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_log_fve_a2_a2mpi_a2mk_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (log(x[i,2] / p[1] ^ 2) + sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3] + p[12] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+12] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_fve_a2(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+9] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_fve_a2_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_fve_a2_a3_a4(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,1] ^ 0.5 + p[11] * x[i,1]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_fve_a2_a2mpi(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+10] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_fve_a2_a2mpi_a2mk(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3]) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_fve_a2_a2mpi_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+11] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end

        function model_plot_ch_fve_a2_a2mpi_a2mk_a3(x,p)
            mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	        nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            return [p[2] + p[3] * x[i,2] / p[1] ^ 2 + p[4] * x[i,3] / p[1] ^ 2 + p[5] * x[i,2] ^ 2 / p[1] ^ 4 + p[6] * x[i,3] ^ 2 / p[1] ^ 4 + p[7] * x[i,2] * x[i,3] / p[1] ^ 4 + p[8] * x[i,2] ^ 4 / p[1] ^ 4 * (sum(mm ./ sqrt.(nn) ./ x[i,4] .* besselk.(1, x[i,4]))) + (p[9] + p[10] * x[i,2] + p[11] * x[i,3] + p[12] * x[i,1] ^ 0.5) * x[i,1] for i in 1:length(x[:,1])]
            #g = [p[i+12] for i in 1:3*div(length(x[:,1]),4)]
            #return [f;g]
        end
    ##

    models = [model_ch_a2, model_ch_a2_a3, model_ch_a2_a2mpi, model_ch_a2_a2mpi_a2mk, model_ch_a2_a2mpi_a3, model_ch_a2_a2mpi_a2mk_a3,
              model_ch_log_a2, model_ch_log_a2_a3, model_ch_log_a2_a2mpi, model_ch_log_a2_a2mpi_a2mk, model_ch_log_a2_a2mpi_a3, model_ch_log_a2_a2mpi_a2mk_a3,
              model_ch_log_fve_a2, model_ch_log_fve_a2_a3, model_ch_log_fve_a2_a2mpi, model_ch_log_fve_a2_a2mpi_a2mk, model_ch_log_fve_a2_a2mpi_a3, model_ch_log_fve_a2_a2mpi_a2mk_a3,
              model_ch_fve_a2, model_ch_fve_a2_a3, model_ch_fve_a2_a2mpi, model_ch_fve_a2_a2mpi_a2mk, model_ch_fve_a2_a2mpi_a3, model_ch_fve_a2_a2mpi_a2mk_a3
    ]
    models_plot = [model_plot_ch_a2, model_plot_ch_a2_a3, model_plot_ch_a2_a2mpi, model_plot_ch_a2_a2mpi_a2mk, model_plot_ch_a2_a2mpi_a3, model_plot_ch_a2_a2mpi_a2mk_a3,
                   model_plot_ch_log_a2, model_plot_ch_log_a2_a3, model_plot_ch_log_a2_a2mpi, model_plot_ch_log_a2_a2mpi_a2mk, model_plot_ch_log_a2_a2mpi_a3, model_plot_ch_log_a2_a2mpi_a2mk_a3,
                   model_plot_ch_log_fve_a2, model_plot_ch_log_fve_a2_a3, model_plot_ch_log_fve_a2_a2mpi, model_plot_ch_log_fve_a2_a2mpi_a2mk, model_plot_ch_log_fve_a2_a2mpi_a3, model_plot_ch_log_fve_a2_a2mpi_a2mk_a3,
                   model_plot_ch_fve_a2, model_plot_ch_fve_a2_a3, model_plot_ch_fve_a2_a2mpi, model_plot_ch_fve_a2_a2mpi_a2mk, model_plot_ch_fve_a2_a2mpi_a3, model_plot_ch_fve_a2_a2mpi_a2mk_a3
    ]
    params = [7,9,9,10,10,11,9,10,10,11,11,12,9,10,10,11,11,12,9,10,10,11,11,12]

    list_all = [
            [ens_340; ens_346; ens_355; ens_370; ens_385],
            [ens_346; ens_355; ens_370; ens_385],
            [2,3,4,6,7,8,9,11,12,13,14,16,17,19],     #nosym
            [1,2,3,4,5,6,7,8,10,11,12,13,15,17,18,19] #mL>4.1
    ]
    for j in 1:length(list_all)
        list = list_all[j]
        xx = value.([x[list,:]; x[list,:]; x[list,:]; x[list,:]])
        yy = [y[list]; x[list,1]; x[list,2]; x[list,3]]

        #i = 1
        for i in 1:length(models)
            uprm, chi2, chi_exp, pval_aux = fit_alg(models[i], xx, yy, params[i]+3*div(length(yy),4), wpm=wpm)
            doff = length(yy) - (params[i]+3*div(length(yy),4))

            function fun(x,p)
                a = [[0.0,p[1],p[2],1000]'; [0.0,p[1],p[2],1000]'; [0.0,p[1],p[2],1000]'; [0.0,p[1],p[2],1000]']
                f = models[i](a,x)[1] .^ 2 ./ p'[:,1]
                g = models[i](a,x)[1] .^ 2 ./ p'[:,2]
                h = [p[2+k] for k in 1:div(length(x),2)]
                return [f;g;h]
            end

            #a = fit_alg(fun, [value.(uprm);value.(uprm)], [ratio1_ph;ratio2_ph;uprm], 2+length(uprm), [value(phi2_ph),value(phi4_ph)])
            #x_ph = [[0.0 a[1][1] a[1][2] 1000]; [0.0 a[1][1] a[1][2] 1000]; [0.0 a[1][1] a[1][2] 1000]; [0.0 a[1][1] a[1][2] 1000]]
            x_ph = [[0.0 phi2_ph phi4_ph 1000]; [0.0 phi2_ph phi4_ph 1000]; [0.0 phi2_ph phi4_ph 1000]; [0.0 phi2_ph phi4_ph 1000]]
            t0mO_ph = models[i](x_ph, uprm)[1]; t0_ph_i = t0mO_ph / mO_ph * hc; uwerr(t0_ph_i); t0_ph_i
            
            push!(t0_ph, t0_ph_i)
            push!(pval_t0, pval_aux)
            push!(dof_t0, doff)
            push!(TIC_t0, chi2 - 2 * chi_exp)

            if i == 1 && j == 1 global uprm_plot = uprm end
        end
    end
    W = exp.(-0.5 * TIC_t0) / sum(exp.(-0.5 * TIC_t0))
    t0_ph_av = sum(t0_ph .* W)
    t02_ph_av = sum(t0_ph .^ 2 .* W)
    syst_t0 = sqrt(t02_ph_av - t0_ph_av ^ 2)
    uwerr(t0_ph_av)
    t0_ph_av
    syst_t0
    #t0_ph_ini = deepcopy(t0_ph)

    uprm = uprm_plot

    function proj_phi2_t0mO()
        Fph = models_plot[i]([1 ./ t0 [phi2_ph for i in 1:length(t0)] phi4 [1000 for i in 1:length(t0)]], uprm)
        Fphi2 = models_plot[i]([1 ./ t0 phi2 phi4 value.(x[:,4])], uprm)
        y_proj = y .- (Fphi2 .- Fph)
        uwerr.(y_proj)

        x_plot = [i for i in 1.08:0.01:1.15]
        #errorbar(value.(x[ens_334,3]), value.(y_proj[ens_334]), err.(y_proj[ens_334]), err.(x[ens_334,3]), fmt="s", label=L"$\beta=3.34$", color="brown")
        errorbar(value.(x[ens_340,3]), value.(y_proj[ens_340]), err.(y_proj[ens_340]), err.(x[ens_340,3]), fmt="s", label=L"$\beta=3.40$", color="purple")
        aux = models_plot[i]([[1 / t0[ens_340[1]] for i in 1:length(x_plot)] [phi2_ph for i in 1:length(x_plot)] x_plot [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="purple", alpha=0.2)
        errorbar(value.(x[ens_346,3]), value.(y_proj[ens_346]), err.(y_proj[ens_346]), err.(x[ens_346,3]), fmt="<", label=L"$\beta=3.46$", color="green")
        aux = models_plot[i]([[1 / t0[ens_346[1]] for i in 1:length(x_plot)] [phi2_ph for i in 1:length(x_plot)] x_plot [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="green", alpha=0.2)
        errorbar(value.(x[ens_355,3]), value.(y_proj[ens_355]), err.(y_proj[ens_355]), err.(x[ens_355,3]), fmt=">", label=L"$\beta=3.55$", color="blue")
        aux = models_plot[i]([[1 / t0[ens_355[1]] for i in 1:length(x_plot)] [phi2_ph for i in 1:length(x_plot)] x_plot [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="blue", alpha=0.2)
        errorbar(value.(x[ens_370,3]), value.(y_proj[ens_370]), err.(y_proj[ens_370]), err.(x[ens_370,3]), fmt="^", label=L"$\beta=3.70$", color="orange")
        aux = models_plot[i]([[1 / t0[ens_370[1]] for i in 1:length(x_plot)] [phi2_ph for i in 1:length(x_plot)] x_plot [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="orange", alpha=0.2)
        errorbar(value.(x[ens_385,3]), value.(y_proj[ens_385]), err.(y_proj[ens_385]), err.(x[ens_385,3]), fmt="o", label=L"$\beta=3.85$", color="red")
        aux = models_plot[i]([[1 / t0[ens_385[1]] for i in 1:length(x_plot)] [phi2_ph for i in 1:length(x_plot)] x_plot [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="red", alpha=0.2)

        #errorbar(value.(x[:,3]), value.(y), err.(y), err.(x[:,3]), fmt="+", label="raw")

        legend()
        xlabel(L"$\phi_4$")
        ylabel(L"$\sqrt{t_0}m_{\Omega}$")
        tight_layout()
    end

    function proj_phi4_t0mO()
        Fph = models_plot[i]([1 ./ t0 phi2 [phi4_ph for i in 1:length(t0)] [1000 for i in 1:length(t0)]], uprm)
        Fphi2 = models_plot[i]([1 ./ t0 phi2 phi4 value.(x[:,4])], uprm)
        y_proj = y .- (Fphi2 .- Fph)
        uwerr.(y_proj)

        x_plot = [i for i in 0.07:0.01:1.1]
        #errorbar(value.(x[ens_334,2]), value.(y_proj[ens_334]), err.(y_proj[ens_334]), err.(x[ens_334,2]), fmt="s", label=L"$\beta=3.34$", color="brown")
        errorbar(value.(x[ens_340,2]), value.(y_proj[ens_340]), err.(y_proj[ens_340]), err.(x[ens_340,2]), fmt="s", label=L"$\beta=3.40$", color="purple")
        aux = models_plot[i]([[1 / t0[ens_340[1]] for i in 1:length(x_plot)] x_plot [phi4_ph for i in 1:length(x_plot)] [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="purple", alpha=0.2)
        errorbar(value.(x[ens_346,2]), value.(y_proj[ens_346]), err.(y_proj[ens_346]), err.(x[ens_346,2]), fmt="<", label=L"$\beta=3.46$", color="green")
        aux = models_plot[i]([[1 / t0[ens_346[1]] for i in 1:length(x_plot)] x_plot [phi4_ph for i in 1:length(x_plot)] [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="green", alpha=0.2)
        errorbar(value.(x[ens_355,2]), value.(y_proj[ens_355]), err.(y_proj[ens_355]), err.(x[ens_355,2]), fmt=">", label=L"$\beta=3.55$", color="blue")
        aux = models_plot[i]([[1 / t0[ens_355[1]] for i in 1:length(x_plot)] x_plot [phi4_ph for i in 1:length(x_plot)] [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="blue", alpha=0.2)
        errorbar(value.(x[ens_370,2]), value.(y_proj[ens_370]), err.(y_proj[ens_370]), err.(x[ens_370,2]), fmt="^", label=L"$\beta=3.70$", color="orange")
        aux = models_plot[i]([[1 / t0[ens_370[1]] for i in 1:length(x_plot)] x_plot [phi4_ph for i in 1:length(x_plot)] [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="orange", alpha=0.2)
        errorbar(value.(x[ens_385,2]), value.(y_proj[ens_385]), err.(y_proj[ens_385]), err.(x[ens_385,2]), fmt="o", label=L"$\beta=3.85$", color="red")
        aux = models_plot[i]([[1 / t0[ens_385[1]] for i in 1:length(x_plot)] x_plot [phi4_ph for i in 1:length(x_plot)] [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux); v = value.(aux); e = err.(aux); fill_between(x_plot, v-e, v+e, color="red", alpha=0.2)
        #errorbar(value.(x[:,2]), value.(y), err.(y), err.(x[:,2]), fmt="+", label="raw")
        legend()
        xlabel(L"$\phi_2$")
        ylabel(L"$\sqrt{t_0}m_{\Omega}$")
        tight_layout()
    end

    function proj_phi24_t0mO()
        Fph = models_plot[i]([1 ./ t0 [phi2_ph for i in 1:length(t0)] [phi4_ph for i in 1:length(t0)] [1000 for i in 1:length(t0)]], uprm)
        Fphi2 = models_plot[i]([1 ./ t0 phi2 phi4 value.(x[:,4])], uprm)
        y_proj = y .- (Fphi2 .- Fph)
        uwerr.(y_proj)

        #errorbar(value.(x[ens_334,1]), value.(y_proj[ens_334]), err.(y_proj[ens_334]), err.(x[ens_334,1]), fmt="s", label=L"$\beta=3.34$", color="brown")
        errorbar(value.(x[ens_340,1]), value.(y_proj[ens_340]), err.(y_proj[ens_340]), err.(x[ens_340,1]), fmt="s", label=L"$\beta=3.40$", color="purple")
        errorbar(value.(x[ens_346,1]), value.(y_proj[ens_346]), err.(y_proj[ens_346]), err.(x[ens_346,1]), fmt="<", label=L"$\beta=3.46$", color="green")
        errorbar(value.(x[ens_355,1]), value.(y_proj[ens_355]), err.(y_proj[ens_355]), err.(x[ens_355,1]), fmt=">", label=L"$\beta=3.55$", color="blue")
        errorbar(value.(x[ens_370,1]), value.(y_proj[ens_370]), err.(y_proj[ens_370]), err.(x[ens_370,1]), fmt="^", label=L"$\beta=3.70$", color="orange")
        errorbar(value.(x[ens_385,1]), value.(y_proj[ens_385]), err.(y_proj[ens_385]), err.(x[ens_385,1]), fmt="o", label=L"$\beta=3.85$", color="red")
        #errorbar(value.(x[:,1]), value.(y), err.(y), err.(x[:,1]), fmt="+", label="raw")
        legend()
        x_plot = [i for i in 0.0:0.01:0.36]
        aux = models_plot[i]([x_plot [phi2_ph for i in 1:length(x_plot)] [phi4_ph for i in 1:length(x_plot)] [1000 for i in 1:length(x_plot)]], uprm) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot, v-e, v+e, color="gray", alpha=0.5)
        legend()
        xlabel(L"$a^2/t_0$")
        ylabel(L"$\sqrt{t_0}m_{\Omega}$")
        tight_layout()
    end

    mods = [
        "[no log][no fve][a2]",
        "[no log][no fve][a2, a3]",
        "[no log][no fve][a2, a2phi2]",
        "[no log][no fve][a2, a2phi2, a2phi4]",
        "[no log][no fve][a2, a2phi2, a3]",
        "[no log][no fve][a2, a2phi2, a2phi4, a3]",
        "[log][no fve][a2]",
        "[log][no fve][a2, a3]",
        "[log][no fve][a2, a2phi2]",
        "[log][no fve][a2, a2phi2, a2phi4]",
        "[log][no fve][a2, a2phi2, a3]",
        "[log][no fve][a2, a2phi2, a2phi4, a3]",
        "[log][fve][a2]",
        "[log][fve][a2, a3]",
        "[log][fve][a2, a2phi2]",
        "[log][fve][a2, a2phi2, a2phi4]",
        "[log][fve][a2, a2phi2, a3]",
        "[log][fve][a2, a2phi2, a2phi4, a3]",
        "[no log][fve][a2]",
        "[no log][fve][a2, a3]",
        "[no log][fve][a2, a2phi2]",
        "[no log][fve][a2, a2phi2, a2phi4]",
        "[no log][fve][a2, a2phi2, a3]",
        "[no log][fve][a2, a2phi2, a2phi4, a3]",

        "[no log][no fve][a2][b3.40]",
        "[no log][no fve][a2, a3][b3.40]",
        "[no log][no fve][a2, a2phi2][b3.40]",
        "[no log][no fve][a2, a2phi2, a2phi4][b3.40]",
        "[no log][no fve][a2, a2phi2, a3][b3.40]",
        "[no log][no fve][a2, a2phi2, a2phi4, a3][b3.40]",
        "[log][no fve][a2][b3.40]",
        "[log][no fve][a2, a3][b3.40]",
        "[log][no fve][a2, a2phi2][b3.40]",
        "[log][no fve][a2, a2phi2, a2phi4][b3.40]",
        "[log][no fve][a2, a2phi2, a3][b3.40]",
        "[log][no fve][a2, a2phi2, a2phi4, a3][b3.40]",
        "[log][fve][a2][b3.40]",
        "[log][fve][a2, a3][b3.40]",
        "[log][fve][a2, a2phi2][b3.40]",
        "[log][fve][a2, a2phi2, a2phi4][b3.40]",
        "[log][fve][a2, a2phi2, a3][b3.40]",
        "[log][fve][a2, a2phi2, a2phi4, a3][b3.40]",
        "[no log][fve][a2][b3.40]",
        "[no log][fve][a2, a3][b3.40]",
        "[no log][fve][a2, a2phi2][b3.40]",
        "[no log][fve][a2, a2phi2, a2phi4][b3.40]",
        "[no log][fve][a2, a2phi2, a3][b3.40]",
        "[no log][fve][a2, a2phi2, a2phi4, a3][b3.40]",

        "[no log][no fve][a2][sym]",
        "[no log][no fve][a2, a3][sym]",
        "[no log][no fve][a2, a2phi2][sym]",
        "[no log][no fve][a2, a2phi2, a2phi4][sym]",
        "[no log][no fve][a2, a2phi2, a3][sym]",
        "[no log][no fve][a2, a2phi2, a2phi4, a3][sym]",
        "[log][no fve][a2][sym]",
        "[log][no fve][a2, a3][sym]",
        "[log][no fve][a2, a2phi2][sym]",
        "[log][no fve][a2, a2phi2, a2phi4][sym]",
        "[log][no fve][a2, a2phi2, a3][sym]",
        "[log][no fve][a2, a2phi2, a2phi4, a3][sym]",
        "[log][fve][a2][sym]",
        "[log][fve][a2, a3][sym]",
        "[log][fve][a2, a2phi2][sym]",
        "[log][fve][a2, a2phi2, a2phi4][sym]",
        "[log][fve][a2, a2phi2, a3][sym]",
        "[log][fve][a2, a2phi2, a2phi4, a3][sym]",
        "[no log][fve][a2][sym]",
        "[no log][fve][a2, a3][sym]",
        "[no log][fve][a2, a2phi2][sym]",
        "[no log][fve][a2, a2phi2, a2phi4][sym]",
        "[no log][fve][a2, a2phi2, a3][sym]",
        "[no log][fve][a2, a2phi2, a2phi4, a3][sym]",

        "[no log][no fve][a2][4.1]",
        "[no log][no fve][a2, a3][4.1]",
        "[no log][no fve][a2, a2phi2][4.1]",
        "[no log][no fve][a2, a2phi2, a2phi4][4.1]",
        "[no log][no fve][a2, a2phi2, a3][4.1]",
        "[no log][no fve][a2, a2phi2, a2phi4, a3][4.1]",
        "[log][no fve][a2][4.1]",
        "[log][no fve][a2, a3][4.1]",
        "[log][no fve][a2, a2phi2][4.1]",
        "[log][no fve][a2, a2phi2, a2phi4][4.1]",
        "[log][no fve][a2, a2phi2, a3][4.1]",
        "[log][no fve][a2, a2phi2, a2phi4, a3][4.1]",
        "[log][fve][a2][4.1]",
        "[log][fve][a2, a3][4.1]",
        "[log][fve][a2, a2phi2][4.1]",
        "[log][fve][a2, a2phi2, a2phi4][4.1]",
        "[log][fve][a2, a2phi2, a3][4.1]",
        "[log][fve][a2, a2phi2, a2phi4, a3][4.1]",
        "[no log][fve][a2][4.1]",
        "[no log][fve][a2, a3][4.1]",
        "[no log][fve][a2, a2phi2][4.1]",
        "[no log][fve][a2, a2phi2, a2phi4][4.1]",
        "[no log][fve][a2, a2phi2, a3][4.1]",
        "[no log][fve][a2, a2phi2, a2phi4, a3][4.1]"
    ]

    fig = figure("model av")
    subplot(211)
    ax = gca() # Get the handle of the current axis
    ylabel(L"$\sqrt{t_0}\; [fm]$")
    xx = collect(1:length(W))
    errorbar(xx, value.(t0_ph), err.(t0_ph), fmt="x")
    v = [value(t0_ph_av) for i in 1:length(t0_ph)]
    e = [sqrt(err(t0_ph_av) ^ 2 + value(syst_t0) ^ 2) for i in 1:length(t0_ph)]
    fill_between(xx, v-e, v+e, color="deepskyblue", alpha=0.3)
    #ax[:set_ylim]([v[1]-e[1], v[1]+e[1]])
    ax[:set_ylim]([0.13, 0.15])
    plt.xticks(xx, collect(1:length(xx)))
    setp(ax.get_xticklabels(),visible=false)
    ax[:set_xlim]([0, length(xx)+1])
    subplot(212)
    ax = gca() # Get the handle of the current axis
    bar(xx, pval_t0, label="pvalue")
    #plot(xx, dof_t0, label="dof")
    bar(xx, W, label="weight")
    for i in div(length(xx), length(list_all)):div(length(xx), length(list_all)):length(xx)
        vlines(i,0,1,ls="--",color="red")
    end
    ax[:set_xlim]([0, length(xx)+1])
    #plt.xticks(xx, collect(1:length(xx)))
    plt.xticks(xx, mods)
    xticks(rotation=90)
    #ax[:set_ylim]([0.0, maximum(W)])
    legend()
    tight_layout()

#===========================================#














using ADerrors, PyPlot
err = ADerrors.err

BKS = uwreal([0.1467,0.00158],1); uwerr(BKS)
Bali = uwreal([0.4098 / sqrt(8), 0.0025 / sqrt(8)],1); uwerr(Bali)
FLAG21 = uwreal([0.14464,0.00087],1); uwerr(FLAG21)
Strass = uwreal([0.1441,0.0013],1); uwerr(Strass)
UK15 = uwreal([0.1511,0.0022],1); uwerr(UK15)
UK14 = uwreal([0.1439,0.0008],1); uwerr(UK14)
BMW12 = uwreal([0.1465,0.0025],1); uwerr(BMW12)

sqrt_t0_st = uwreal([0.1433,0.0010],""); uwerr(sqrt_t0_st)
sqrt_t0_tm = uwreal([0.1442,0.0011],""); uwerr(sqrt_t0_tm)
sqrt_t0_comb = uwreal([0.1438,0.0008],""); uwerr(sqrt_t0_comb)

sqrt_t0_st_16 = uwreal([0.1437,0.0011],""); uwerr(sqrt_t0_st_16)
sqrt_t0_tm_16 = uwreal([0.1450,0.0011],""); uwerr(sqrt_t0_tm_16)
sqrt_t0_comb_16 = uwreal([0.1436,0.0009],""); uwerr(sqrt_t0_comb_16)

sqrt_t0_st_old = uwreal([0.1444,0.0009],""); uwerr(sqrt_t0_st_old)
sqrt_t0_tm_old = uwreal([0.1445,0.0013],""); uwerr(sqrt_t0_tm_old)
sqrt_t0_comb_old = uwreal([0.1448,0.0007],""); uwerr(sqrt_t0_comb_old)

sqrt_t0_st_old_16 = uwreal([0.1449,0.0010],""); uwerr(sqrt_t0_st_old_16)
sqrt_t0_tm_old_16 = uwreal([0.1449,0.0013],""); uwerr(sqrt_t0_tm_old_16)
sqrt_t0_comb_old_16 = uwreal([0.1453,0.0008],""); uwerr(sqrt_t0_comb_old_16)

fig = figure("",figsize=(5,5))
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 20

#errorbar(value.([sqrt_t0_st_SU2_fpi]), 1, [0], err.([sqrt_t0_st_SU2_fpi]), fmt="s", color="deeppink")
#errorbar(value.([sqrt_t0_st_SU2_fpik]), 3, [0], err.([sqrt_t0_st_SU2_fpik]), fmt="s", color="deeppink")

#errorbar(value.([sqrt_t0_st_old]), 7, [0], err.([sqrt_t0_st_old]), fmt="s", color="green")

errorbar(value.([BMW12]), 1, [0], err.([BMW12]), fmt="s", color="blue")
errorbar(value.([UK14]), 4, [0], err.([UK14]), fmt="s", color="blue")
errorbar(value.([UK15]), 7, [0], err.([UK15]), fmt="s", color="blue")
errorbar(value.([BKS]), 10, [0], err.([BKS]), fmt="s", color="blue")
errorbar(value.([FLAG21]), 13, [0], err.([FLAG21]), fmt="s", color="black")
errorbar(value.([Bali]), 16, [0], err.([Bali]), fmt="s", color="blue")
errorbar(value.([Strass]), 19, [0], err.([Strass]), fmt="s", color="blue")
errorbar(value.([sqrt_t0_st]), 22, [0], err.([sqrt_t0_st]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_tm]), 25, [0], err.([sqrt_t0_tm]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_comb]), 28, [0], err.([sqrt_t0_comb]), fmt="s", color="deeppink")

xlabel(L"$\sqrt{t_0}\;\rm{[fm]}$")
yticks([1,4,7,10,13,16,19,22,25,28], 
       [L"BMW '12 $[m_{\Omega}]$",
        L"RBC/UKQCD '14 $[m_{\Omega}]$",
        L"QCDSF/UKQCD '15 $[m_O,\;M_V]$",
        L"Bruno et al. '16 $[f_{\pi K}]$",
        "FLAG '21", 
        L"Bali et al. '22 $[m_{\Xi}]$",
        L"Strassberger '23 $[f_{\pi K}]$", 
        L"This work, Wilson $[f_{\pi K}]$", 
        L"This work, Wtm $[f_{\pi K}]$",  
        L"This work, Combined $[f_{\pi K}]$"])
x_plot = [FLAG21 for i in 1:1:28]
v = [i for i in 1:1:28]

fill_betweenx(v, value(FLAG21)-err(FLAG21), value(FLAG21)+err(FLAG21), color="lightsteelblue", alpha=0.40)
ax = gca()
ax[:set_xlim]([0.14, 0.154])
tight_layout()
savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0_compar.pdf")

#===============================================================================================================#

fig = figure("FLAG",figsize=(5,5))
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 20

errorbar(value.([sqrt_t0_st_old_16]), 1, [0], err.([sqrt_t0_st_old_16]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_st_old]), 2, [0], err.([sqrt_t0_st_old]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_st]), 3, [0], err.([sqrt_t0_st]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_tm_old_16]), 7, [0], err.([sqrt_t0_tm_old_16]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_tm_16]), 8, [0], err.([sqrt_t0_tm_16]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_tm]), 9, [0], err.([sqrt_t0_tm]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_comb_old_16]), 13, [0], err.([sqrt_t0_comb_old_16]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_comb_16]), 14, [0], err.([sqrt_t0_comb_16]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_comb]), 15, [0], err.([sqrt_t0_comb]), fmt="s", color="deeppink")
errorbar(value.([BKS]), 19, [0], err.([BKS]), fmt="s", color="blue")
errorbar(value.([FLAG21]), 21, [0], err.([FLAG21]), fmt="s", color="black")

xlabel(L"$\sqrt{t_0}\;\rm{[fm]}$")
yticks([1,7,13,19,21], 
       [L"This work, Wilson $[f_{\pi K}]$", 
        L"This work, Wtm $[f_{\pi K}]$",  
        L"This work, Combined $[f_{\pi K}]$",
        L"Bruno et al. '16 $[f_{\pi K}]$",
        "FLAG '21"])
x_plot = [FLAG21 for i in 1:1:22]
v = [i for i in 1:1:22]

fill_betweenx(v, value(FLAG21)-err(FLAG21), value(FLAG21)+err(FLAG21), color="lightsteelblue", alpha=0.40)
ax = gca()
ax[:set_xlim]([0.14, 0.154])
tight_layout()
savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0_compar_FLAG.pdf")
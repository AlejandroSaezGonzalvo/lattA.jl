using ADerrors, PyPlot
err = ADerrors.err

BKS = uwreal([0.1467,0.00158],""); uwerr(BKS)
Moh = uwreal([0.1448,0.0003],""); uwerr(Moh)
Bali = uwreal([0.4098 / sqrt(8), 0.0025 / sqrt(8)],""); uwerr(Bali)
FLAG21 = uwreal([0.14464,0.00087],""); uwerr(FLAG21)
Strass = uwreal([0.1441,0.0013],""); uwerr(Strass)
UK15 = uwreal([0.1511,0.0022],""); uwerr(UK15)
UK14 = uwreal([0.1439,0.0008],""); uwerr(UK14)
BMW12 = uwreal([0.1465,0.0025],""); uwerr(BMW12)

sqrt_t0_st = uwreal([0.1436,0.0009],""); uwerr(sqrt_t0_st)
sqrt_t0_tm = uwreal([0.1444,0.0010],""); uwerr(sqrt_t0_tm)
sqrt_t0_comb = uwreal([0.1441,0.0007],""); uwerr(sqrt_t0_comb)

sqrt_t0_st_AIC = uwreal([0.1433,0.0012],""); uwerr(sqrt_t0_st_AIC) ## this is w/o penal for cuts and with c-switch
sqrt_t0_tm_AIC = uwreal([0.1444,0.0009],""); uwerr(sqrt_t0_tm_AIC)
sqrt_t0_comb_AIC = uwreal([0.1443,0.0007],""); uwerr(sqrt_t0_comb_AIC)

sqrt_t0_st_AIC_sub = uwreal([0.1438,0.0009],""); uwerr(sqrt_t0_st_AIC_sub)
sqrt_t0_tm_AIC_sub = uwreal([0.1446,0.0008],""); uwerr(sqrt_t0_tm_AIC_sub)
sqrt_t0_comb_AIC_sub = uwreal([0.1443,0.0007],""); uwerr(sqrt_t0_comb_AIC_sub)

sqrt_t0_st_AIC_noc = uwreal([0.1432,0.0013],""); uwerr(sqrt_t0_st_AIC) ## this is w/o penal for cuts and w/o c-switch
sqrt_t0_tm_AIC_noc = uwreal([0.1442,0.0009],""); uwerr(sqrt_t0_tm_AIC)
sqrt_t0_comb_AIC_noc = uwreal([0.1442,0.0006],""); uwerr(sqrt_t0_comb_AIC)

sqrt_t0_st_AIC_sub_noc = uwreal([0.1440,0.0009],""); uwerr(sqrt_t0_st_AIC_sub) 
sqrt_t0_tm_AIC_sub_noc = uwreal([0.1443,0.0010],""); uwerr(sqrt_t0_tm_AIC_sub)
sqrt_t0_comb_AIC_sub_noc = uwreal([0.1443,0.0007],""); uwerr(sqrt_t0_comb_AIC_sub)

sqrt_t0_st_noc = uwreal([0.1440,0.0008],""); uwerr(sqrt_t0_st_noc)
sqrt_t0_tm_noc = uwreal([0.1442,0.0009],""); uwerr(sqrt_t0_tm_noc)
sqrt_t0_comb_noc = uwreal([0.1443,0.0006],""); uwerr(sqrt_t0_comb_noc)

sqrt_t0_st_cinf = uwreal([0.1434,0.0012],""); uwerr(sqrt_t0_st_cinf)
sqrt_t0_tm_cinf = uwreal([0.1444,0.0013],""); uwerr(sqrt_t0_tm_cinf)
sqrt_t0_comb_cinf = uwreal([0.1439,0.0009],""); uwerr(sqrt_t0_comb_cinf)

sqrt_t0_st_fpi = uwreal([0.1433,0.0012],""); uwerr(sqrt_t0_st_fpi)

sqrt_t0_sym = uwreal([0.1438,0.0007],""); uwerr(sqrt_t0_sym)
sqrt_t0_star = uwreal([0.1438,0.0008],""); uwerr(sqrt_t0_star)

fK_pred = uwreal([156.7,1.3],""); uwerr(fK_pred)
fK_FLAG21 = uwreal([157.2,0.5],""); uwerr(fK_FLAG21)
fK_FLAG16 = uwreal([156.2,0.7],""); uwerr(fK_FLAG16)

sqrt_t0_st_16 = uwreal([0.1443,0.0011],""); uwerr(sqrt_t0_st_16)
sqrt_t0_tm_16 = uwreal([0.1450,0.0012],""); uwerr(sqrt_t0_tm_16)
sqrt_t0_comb_16 = uwreal([0.1449,0.0008],""); uwerr(sqrt_t0_comb_16)

sqrt_t0_st_AP = uwreal([0.1437,0.0010],""); uwerr(sqrt_t0_st_AP)
sqrt_t0_tm_AP = uwreal([0.1446,0.0010],""); uwerr(sqrt_t0_tm_AP)
sqrt_t0_comb_AP = uwreal([0.1443,0.0007],""); uwerr(sqrt_t0_comb_AP)

#=
sqrt_t0_st_old = uwreal([0.1444,0.0009],""); uwerr(sqrt_t0_st_old)
sqrt_t0_tm_old = uwreal([0.1445,0.0013],""); uwerr(sqrt_t0_tm_old)
sqrt_t0_comb_old = uwreal([0.1448,0.0007],""); uwerr(sqrt_t0_comb_old)

sqrt_t0_st_old_16 = uwreal([0.1449,0.0010],""); uwerr(sqrt_t0_st_old_16)
sqrt_t0_tm_old_16 = uwreal([0.1449,0.0013],""); uwerr(sqrt_t0_tm_old_16)
sqrt_t0_comb_old_16 = uwreal([0.1453,0.0008],""); uwerr(sqrt_t0_comb_old_16)
=#

#=====================================================================================#

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
errorbar(value.([Moh]), 22, [0], err.([Moh]), fmt="s", color="blue")
errorbar(value.([sqrt_t0_st]), 25, [0], err.([sqrt_t0_st]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_tm]), 28, [0], err.([sqrt_t0_tm]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_comb]), 31, [0], err.([sqrt_t0_comb]), fmt="s", color="deeppink")

xlabel(L"$\sqrt{t_0}\;\rm{[fm]}$")
yticks([1,4,7,10,13,16,19,22,25,28,31], 
       [L"BMW '12 $[m_{\Omega}]$",
        L"RBC/UKQCD '14 $[m_{\Omega}]$",
        L"QCDSF/UKQCD '15 $[m_O,\;M_V]$",
        L"Bruno et al. '16 $[f_{\pi K}]$",
        "FLAG '21", 
        L"Bali et al. '22 $[m_{\Xi}]$",
        L"Strassberger '23 $[f_{\pi K}]$",
        L"Hudspith et al. '24 $[m_{\Omega}]$", 
        L"This work, Wilson $[f_{\pi K}]$", 
        L"This work, Wtm $[f_{\pi K}]$",  
        L"This work, Combined $[f_{\pi K}]$"])
x_plot = [FLAG21 for i in 1:1:31]
v = [i for i in 1:1:31]

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

#===============================================================================================================#

fig = figure("swicths",figsize=(5,5))
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 20

errorbar(value.([sqrt_t0_st_noc]), 1, [0], err.([sqrt_t0_st_noc]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_st]), 2, [0], err.([sqrt_t0_st]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_st_cinf]), 3, [0], err.([sqrt_t0_st_cinf]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_tm_noc]), 7, [0], err.([sqrt_t0_tm_noc]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_tm]), 8, [0], err.([sqrt_t0_tm]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_tm_cinf]), 9, [0], err.([sqrt_t0_tm_cinf]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_comb_noc]), 13, [0], err.([sqrt_t0_comb_noc]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_comb]), 14, [0], err.([sqrt_t0_comb]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_comb_cinf]), 15, [0], err.([sqrt_t0_comb_cinf]), fmt="s", color="deeppink")
errorbar(value.([FLAG21]), 19, [0], err.([FLAG21]), fmt="s", color="black")

xlabel(L"$\sqrt{t_0}\;\rm{[fm]}$")
yticks([1,7,13,19], 
       [L"This work, Wilson $[f_{\pi K}]$", 
        L"This work, Wtm $[f_{\pi K}]$",  
        L"This work, Combined $[f_{\pi K}]$",
        "FLAG '21"])
x_plot = [FLAG21 for i in 1:1:20]
v = [i for i in 1:1:20]

fill_betweenx(v, value(FLAG21)-err(FLAG21), value(FLAG21)+err(FLAG21), color="lightsteelblue", alpha=0.40)
xlim(0.14,0.147)
tight_layout()
savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0_compar_switch.pdf")

#===============================================================================================================#


fig = figure("FLAG",figsize=(5,5))
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 20

errorbar(value.([sqrt_t0_st_AIC]), 1, [0], err.([sqrt_t0_st_AIC]), fmt="<", color="deeppink")
errorbar(value.([sqrt_t0_st_AIC_sub]), 2, [0], err.([sqrt_t0_st_AIC_sub]), fmt=">", color="deeppink")
errorbar(value.([sqrt_t0_st]), 3, [0], err.([sqrt_t0_st]), fmt="x", color="deeppink")
errorbar(value.([sqrt_t0_tm_AIC]), 7, [0], err.([sqrt_t0_tm_AIC]), fmt="<", color="deeppink")
errorbar(value.([sqrt_t0_tm_AIC_sub]), 8, [0], err.([sqrt_t0_tm_AIC_sub]), fmt=">", color="deeppink")
errorbar(value.([sqrt_t0_tm]), 9, [0], err.([sqrt_t0_tm]), fmt="x", color="deeppink")
errorbar(value.([sqrt_t0_comb_AIC]), 13, [0], err.([sqrt_t0_comb_AIC]), fmt="<", color="deeppink")
errorbar(value.([sqrt_t0_comb_AIC_sub]), 14, [0], err.([sqrt_t0_comb_AIC_sub]), fmt=">", color="deeppink")
errorbar(value.([sqrt_t0_comb]), 15, [0], err.([sqrt_t0_comb]), fmt="x", color="deeppink")
#errorbar(value.([BKS]), 19, [0], err.([BKS]), fmt="s", color="blue")
errorbar(value.([FLAG21]), 19, [0], err.([FLAG21]), fmt="s", color="black")

xlabel(L"$\sqrt{t_0}\;\rm{[fm]}$")
yticks([1,7,13,19], 
       [L"This work, Wilson $[f_{\pi K}]$", 
        L"This work, Wtm $[f_{\pi K}]$",  
        L"This work, Combined $[f_{\pi K}]$",
        "FLAG '21"])
x_plot = [FLAG21 for i in 1:1:22]
v = [i for i in 1:1:22]

fill_betweenx(v, value(FLAG21)-err(FLAG21), value(FLAG21)+err(FLAG21), color="lightsteelblue", alpha=0.40)
xlim(0.141,0.147)
tight_layout()


savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0_compar_IC.pdf")

#===============================================================================================================#

fig = figure("FLAG",figsize=(5,5))
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 20

errorbar(value.([sqrt_t0_st_16]), 1, [0], err.([sqrt_t0_st_16]), fmt="<", color="deeppink")
errorbar(value.([sqrt_t0_st_AP]), 2, [0], err.([sqrt_t0_st_AP]), fmt=">", color="deeppink")
errorbar(value.([sqrt_t0_st]), 3, [0], err.([sqrt_t0_st]), fmt="x", color="deeppink")
errorbar(value.([sqrt_t0_tm_16]), 7, [0], err.([sqrt_t0_tm_16]), fmt="<", color="deeppink")
errorbar(value.([sqrt_t0_tm_AP]), 8, [0], err.([sqrt_t0_tm_AP]), fmt=">", color="deeppink")
errorbar(value.([sqrt_t0_tm]), 9, [0], err.([sqrt_t0_tm]), fmt="x", color="deeppink")
errorbar(value.([sqrt_t0_comb_16]), 13, [0], err.([sqrt_t0_comb_16]), fmt="<", color="deeppink")
errorbar(value.([sqrt_t0_comb_AP]), 14, [0], err.([sqrt_t0_comb_AP]), fmt=">", color="deeppink")
errorbar(value.([sqrt_t0_comb]), 15, [0], err.([sqrt_t0_comb]), fmt="x", color="deeppink")
#errorbar(value.([BKS]), 19, [0], err.([BKS]), fmt="s", color="blue")
errorbar(value.([FLAG21]), 19, [0], err.([FLAG21]), fmt="s", color="black")

xlabel(L"$\sqrt{t_0}\;\rm{[fm]}$")
yticks([1,7,13,19], 
       [L"This work, Wilson $[f_{\pi K}]$", 
        L"This work, Wtm $[f_{\pi K}]$",  
        L"This work, Combined $[f_{\pi K}]$",
        "FLAG '21"])
x_plot = [FLAG21 for i in 1:1:22]
v = [i for i in 1:1:22]

fill_betweenx(v, value(FLAG21)-err(FLAG21), value(FLAG21)+err(FLAG21), color="lightsteelblue", alpha=0.40)
ax = gca()
ax[:set_xlim]([0.14, 0.154])
tight_layout()
savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0_compar_FLAG_AP.pdf")


#===============================================================================================================#

fig = figure("FLAG",figsize=(5,5))
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 20

errorbar(value.([sqrt_t0_star]), 1, [0], err.([sqrt_t0_star]), fmt="<", color="deeppink")
errorbar(value.([sqrt_t0_sym]), 7, [0], err.([sqrt_t0_sym]), fmt=">", color="deeppink")
errorbar(value.([sqrt_t0_comb]), 13, [0], err.([sqrt_t0_comb]), fmt="x", color="deeppink")
errorbar(value.([FLAG21]), 19, [0], err.([FLAG21]), fmt="s", color="black")

xlabel(L"$\sqrt{t_0}\;\rm{[fm]}$")
yticks([1,7,13,19], 
       ["This work, * point", 
        "This work, sym. point",  
        "This work, ph. point",
        "FLAG '21"])
x_plot = [FLAG21 for i in 0:1:20]
v = [i for i in 0:1:20]

fill_betweenx(v, value(FLAG21)-err(FLAG21), value(FLAG21)+err(FLAG21), color="lightsteelblue", alpha=0.40)
ax = gca()
ax[:set_xlim]([0.142, 0.147])
tight_layout()
savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0_compar_syms.pdf")


#===============================================================================================================#

fig = figure("FLAG",figsize=(5,5))
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 20

errorbar(value.([sqrt_t0_st_fpi]), 1, [0], err.([sqrt_t0_st_fpi]), fmt="<", color="deeppink")
errorbar(value.([sqrt_t0_st]), 7, [0], err.([sqrt_t0_st]), fmt=">", color="deeppink")
errorbar(value.([sqrt_t0_comb]), 13, [0], err.([sqrt_t0_comb]), fmt="x", color="deeppink")
errorbar(value.([FLAG21]), 19, [0], err.([FLAG21]), fmt="s", color="black")

xlabel(L"$\sqrt{t_0}\;\rm{[fm]}$")
yticks([1,7,13,19], 
       [L"This work, Wilson $[f_{\pi}]$", 
        L"This work, Wilson $[f_{\pi K}]$",  
        L"This work, Combined $[f_{\pi K}]$",
        "FLAG '21"])
x_plot = [FLAG21 for i in 0:1:20]
v = [i for i in 0:1:20]

fill_betweenx(v, value(FLAG21)-err(FLAG21), value(FLAG21)+err(FLAG21), color="lightsteelblue", alpha=0.40)
ax = gca()
ax[:set_xlim]([0.142, 0.147])
tight_layout()
savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0_compar_fpi.pdf")


fig = figure("FLAG",figsize=(5,5))
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 20

errorbar(value.([fK_pred]), 1, [0], err.([fK_pred]), fmt="<", color="deeppink")
errorbar(value.([fK_FLAG16]), 7, [0], err.([fK_FLAG16]), fmt="x", color="black")
errorbar(value.([fK_FLAG21]), 13, [0], err.([fK_FLAG21]), fmt="x", color="black")

xlabel(L"$f_K\;\rm{[MeV]}$")
yticks([1,7,13], 
       [L"This work, Wilson, $t_0[f_{\pi}]$", 
        "FLAG '16",
        "FLAG '21"])
x_plot = [fK_FLAG21 for i in 0:1:14]
v = [i for i in 0:1:14]

fill_betweenx(v, value(fK_FLAG21)-err(fK_FLAG21), value(fK_FLAG21)+err(fK_FLAG21), color="lightsteelblue", alpha=0.40)
tight_layout()
savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/fk_compar.pdf")


#===============================================================================================================#


using juobs, BDIO, ADerrors, PyPlot
err = ADerrors.err

## read
        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_1_0.0035_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_1_0035_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_1_0035_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_1_0035_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_0_0_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_0_0_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_0_0_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_0_0_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_10_0_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_10_0_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_10_0_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_10_0_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_10_0.0035_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_10_0035_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_10_0035_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_10_0035_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_10_10_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_10_10_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_10_10_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_10_10_21_c = read_uwreal(fb)







        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_1_0.0035_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_1_0035_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_1_0035_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_1_0035_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_0_0_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_0_0_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_0_0_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_0_0_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_10_0_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_10_0_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_10_0_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_10_0_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_10_0.0035_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_10_0035_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_10_0035_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_10_0035_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_10_10_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_10_10_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_10_10_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_10_10_16_c = read_uwreal(fb)








        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_1_0.0035_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_1_0035_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_1_0035_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_1_0035_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_0_0_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_0_0_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_0_0_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_0_0_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_10_0_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_10_0_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_10_0_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_10_0_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_10_0.0035_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_10_0035_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_10_0035_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_10_0035_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_10_10_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_10_10_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_10_10_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_10_10_21_c = read_uwreal(fb)







        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_1_0.0035_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_1_0035_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_1_0035_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_1_0035_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_0_0_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_0_0_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_0_0_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_0_0_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_10_0_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_10_0_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_10_0_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_10_0_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_10_0.0035_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_10_0035_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_10_0035_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_10_0035_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_10_10_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_10_10_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_10_10_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_10_10_16_c = read_uwreal(fb)
##

t0 = [t0_all_1_0035_21_wtm, t0_all_1_0035_21_wil, t0_all_1_0035_21_c,
      t0_all_0_0_21_wtm, t0_all_0_0_21_wil, t0_all_0_0_21_c,
      t0_all_10_0_21_wtm, t0_all_10_0_21_wil, t0_all_10_0_21_c,
      t0_all_10_0035_21_wtm, t0_all_10_0035_21_wil, t0_all_10_0035_21_c,
      t0_all_10_10_21_wtm, t0_all_10_10_21_wil, t0_all_10_10_21_c,

      t0_all_1_0035_16_wtm, t0_all_1_0035_16_wil, t0_all_1_0035_16_c,
      t0_all_0_0_16_wtm, t0_all_0_0_16_wil, t0_all_0_0_16_c,
      t0_all_10_0_16_wtm, t0_all_10_0_16_wil, t0_all_10_0_16_c,
      t0_all_10_0035_16_wtm, t0_all_10_0035_16_wil, t0_all_10_0035_16_c,
      t0_all_10_10_16_wtm, t0_all_10_10_16_wil, t0_all_10_10_16_c,

      t0_old_1_0035_21_wtm, t0_old_1_0035_21_wil, t0_old_1_0035_21_c,
      t0_old_0_0_21_wtm, t0_old_0_0_21_wil, t0_old_0_0_21_c,
      t0_old_10_0_21_wtm, t0_old_10_0_21_wil, t0_old_10_0_21_c,
      t0_old_10_0035_21_wtm, t0_old_10_0035_21_wil, t0_old_10_0035_21_c,
      t0_old_10_10_21_wtm, t0_old_10_10_21_wil, t0_old_10_10_21_c,

      t0_old_1_0035_16_wtm, t0_old_1_0035_16_wil, t0_old_1_0035_16_c,
      t0_old_0_0_16_wtm, t0_old_0_0_16_wil, t0_old_0_0_16_c,
      t0_old_10_0_16_wtm, t0_old_10_0_16_wil, t0_old_10_0_16_c,
      t0_old_10_0035_16_wtm, t0_old_10_0035_16_wil, t0_old_10_0035_16_c,
      t0_old_10_10_16_wtm, t0_old_10_10_16_wil, t0_old_10_10_16_c
]

sqrt_t0 = sqrt.(t0); uwerr.(sqrt_t0)
sqrt_t0_wtm = sqrt_t0[1:3:end]
sqrt_t0_wil = sqrt_t0[2:3:end]
sqrt_t0_c = sqrt_t0[3:3:end]

BKS = uwreal([0.1467,0.00158],1); uwerr(BKS)
Bali = uwreal([0.4098 / sqrt(8), 0.0025 / sqrt(8)],1); uwerr(Bali)
FLAG21 = uwreal([0.14464,0.00087],1); uwerr(FLAG21)
Strass = uwreal([0.1441,0.0013],1); uwerr(Strass)
UK15 = uwreal([0.1511,0.0022],1); uwerr(UK15)
UK14 = uwreal([0.1439,0.0008],1); uwerr(UK14)
BMW12 = uwreal([0.1465,0.0025],1); uwerr(BMW12)

fig = figure("individual",figsize=(5,5))
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 10

errorbar(value.([BMW12]), 1, [0], err.([BMW12]), fmt="s", color="blue")
errorbar(value.([UK14]), 2, [0], err.([UK14]), fmt="s", color="blue")
errorbar(value.([UK15]), 3, [0], err.([UK15]), fmt="s", color="blue")
errorbar(value.([BKS]), 4, [0], err.([BKS]), fmt="s", color="blue")
errorbar(value.([FLAG21]), 5, [0], err.([FLAG21]), fmt="s", color="black")
errorbar(value.([Bali]), 6, [0], err.([Bali]), fmt="s", color="blue")
errorbar(value.([Strass]), 7, [0], err.([Strass]), fmt="s", color="blue")
for i in 1:length(sqrt_t0_c)
        errorbar(value(sqrt_t0_c[i]), 7+i, 0, err(sqrt_t0_c[i]), fmt="s", color="deeppink")
end
v = [i for i in 1:1:28]
fill_betweenx(v, value(FLAG21)-err(FLAG21), value(FLAG21)+err(FLAG21), color="lightsteelblue", alpha=0.40)

yticks(collect(1:7+length(sqrt_t0_c)), 
       [L"BMW '12 $[m_{\Omega}]$",
        L"RBC/UKQCD '14 $[m_{\Omega}]$",
        L"QCDSF/UKQCD '15 $[m_O,\;M_V]$",
        L"Bruno et al. '16 $[f_{\pi K}]$",
        "FLAG '21", 
        L"Bali et al. '22 $[m_{\Xi}]$",
        L"Strassberger '23 $[f_{\pi K}]$",

        L"all, FLAG21, $c_{\beta}=1,\;c_{\phi_2}=0.0035$",
        L"all, FLAG21, $c_{\beta}=c_{\phi_2}=0$",
        L"all, FLAG21, $c_{\beta}=10,\;c_{\phi_2}=0$",
        L"all, FLAG21, $c_{\beta}=10,\;c_{\phi_2}=0.0035$",
        L"all, FLAG21, $c_{\beta}=10,\;c_{\phi_2}=10$",

        L"all, FLAG16, $c_{\beta}=1,\;c_{\phi_2}=0.0035$",
        L"all, FLAG16, $c_{\beta}=c_{\phi_2}=0$",
        L"all, FLAG16, $c_{\beta}=10,\;c_{\phi_2}=0$",
        L"all, FLAG16, $c_{\beta}=10,\;c_{\phi_2}=0.0035$",
        L"all, FLAG16, $c_{\beta}=10,\;c_{\phi_2}=10$",

        L"old, FLAG21, $c_{\beta}=1,\;c_{\phi_2}=0.0035$",
        L"old, FLAG21, $c_{\beta}=c_{\phi_2}=0$",
        L"old, FLAG21, $c_{\beta}=10,\;c_{\phi_2}=0$",
        L"old, FLAG21, $c_{\beta}=10,\;c_{\phi_2}=0.0035$",
        L"old, FLAG21, $c_{\beta}=10,\;c_{\phi_2}=10$",

        L"old, FLAG16, $c_{\beta}=1,\;c_{\phi_2}=0.0035$",
        L"old, FLAG16, $c_{\beta}=c_{\phi_2}=0$",
        L"old, FLAG16, $c_{\beta}=10,\;c_{\phi_2}=0$",
        L"old, FLAG16, $c_{\beta}=10,\;c_{\phi_2}=0.0035$",
        L"old, FLAG16, $c_{\beta}=10,\;c_{\phi_2}=10$"
])

xlabel(L"$\sqrt{t_0}\;\rm{[fm]}$")
ax = gca()
ax[:set_xlim]([0.14, 0.154])
plt.title("Combined")
tight_layout()

savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0_variations_comb.pdf")

#====================================================================================#

using juobs, BDIO, ADerrors, PyPlot
err = ADerrors.err

plt.ioff()

## read
        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_1_0.0035_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_1_0035_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_1_0035_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_1_0035_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_0_0_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_0_0_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_0_0_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_0_0_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_10_0_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_10_0_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_10_0_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_10_0_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_10_0.0035_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_10_0035_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_10_0035_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_10_0035_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_10_10_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_10_10_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_10_10_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_10_10_21_c = read_uwreal(fb)







        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_1_0.0035_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_1_0035_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_1_0035_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_1_0035_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_0_0_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_0_0_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_0_0_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_0_0_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_10_0_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_10_0_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_10_0_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_10_0_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_10_0.0035_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_10_0035_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_10_0035_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_10_0035_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_allens_10_10_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_all_10_10_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_all_10_10_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_all_10_10_16_c = read_uwreal(fb)








        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_1_0.0035_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_1_0035_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_1_0035_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_1_0035_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_0_0_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_0_0_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_0_0_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_0_0_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_10_0_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_10_0_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_10_0_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_10_0_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_10_0.0035_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_10_0035_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_10_0035_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_10_0035_21_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_10_10_FLAG21.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_10_10_21_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_10_10_21_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_10_10_21_c = read_uwreal(fb)







        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_1_0.0035_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_1_0035_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_1_0035_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_1_0035_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_0_0_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_0_0_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_0_0_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_0_0_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_10_0_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_10_0_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_10_0_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_10_0_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_10_0.0035_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_10_0035_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_10_0035_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_10_0035_16_c = read_uwreal(fb)

        fb = BDIO_open("/home/asaez/cls_ens/results/t0_variations/t0_oldens_10_10_FLAG16.bdio", "r") 
        BDIO_seek!(fb) 
        t0_old_10_10_16_wtm = read_uwreal(fb)
        BDIO_seek!(fb, 2)
        t0_old_10_10_16_wil = read_uwreal(fb)
        BDIO_seek!(fb, 2)  
        t0_old_10_10_16_c = read_uwreal(fb)
##

t0 = [t0_all_1_0035_21_wtm, t0_all_1_0035_21_wil, t0_all_1_0035_21_c,
      t0_all_0_0_21_wtm, t0_all_0_0_21_wil, t0_all_0_0_21_c,
      t0_all_10_0_21_wtm, t0_all_10_0_21_wil, t0_all_10_0_21_c,
      t0_all_10_0035_21_wtm, t0_all_10_0035_21_wil, t0_all_10_0035_21_c,
      t0_all_10_10_21_wtm, t0_all_10_10_21_wil, t0_all_10_10_21_c,

      t0_all_1_0035_16_wtm, t0_all_1_0035_16_wil, t0_all_1_0035_16_c,
      t0_all_0_0_16_wtm, t0_all_0_0_16_wil, t0_all_0_0_16_c,
      t0_all_10_0_16_wtm, t0_all_10_0_16_wil, t0_all_10_0_16_c,
      t0_all_10_0035_16_wtm, t0_all_10_0035_16_wil, t0_all_10_0035_16_c,
      t0_all_10_10_16_wtm, t0_all_10_10_16_wil, t0_all_10_10_16_c,

      t0_old_1_0035_21_wtm, t0_old_1_0035_21_wil, t0_old_1_0035_21_c,
      t0_old_0_0_21_wtm, t0_old_0_0_21_wil, t0_old_0_0_21_c,
      t0_old_10_0_21_wtm, t0_old_10_0_21_wil, t0_old_10_0_21_c,
      t0_old_10_0035_21_wtm, t0_old_10_0035_21_wil, t0_old_10_0035_21_c,
      t0_old_10_10_21_wtm, t0_old_10_10_21_wil, t0_old_10_10_21_c,

      t0_old_1_0035_16_wtm, t0_old_1_0035_16_wil, t0_old_1_0035_16_c,
      t0_old_0_0_16_wtm, t0_old_0_0_16_wil, t0_old_0_0_16_c,
      t0_old_10_0_16_wtm, t0_old_10_0_16_wil, t0_old_10_0_16_c,
      t0_old_10_0035_16_wtm, t0_old_10_0035_16_wil, t0_old_10_0035_16_c,
      t0_old_10_10_16_wtm, t0_old_10_10_16_wil, t0_old_10_10_16_c
]

sqrt_t0 = sqrt.(t0); uwerr.(sqrt_t0)
sqrt_t0_wtm = sqrt_t0[1:3:end]
sqrt_t0_wil = sqrt_t0[2:3:end]
sqrt_t0_c = sqrt_t0[3:3:end]

BKS = uwreal([0.1467,0.00158],1); uwerr(BKS)
Bali = uwreal([0.4098 / sqrt(8), 0.0025 / sqrt(8)],1); uwerr(Bali)
FLAG21 = uwreal([0.14464,0.00087],1); uwerr(FLAG21)
Strass = uwreal([0.1441,0.0013],1); uwerr(Strass)
UK15 = uwreal([0.1511,0.0022],1); uwerr(UK15)
UK14 = uwreal([0.1439,0.0008],1); uwerr(UK14)
BMW12 = uwreal([0.1465,0.0025],1); uwerr(BMW12)

fig = figure("together",figsize=(10,15))
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 10

errorbar(value.([BMW12]), 1, [0], err.([BMW12]), fmt="s", color="blue")
errorbar(value.([UK14]), 2, [0], err.([UK14]), fmt="s", color="blue")
errorbar(value.([UK15]), 3, [0], err.([UK15]), fmt="s", color="blue")
errorbar(value.([BKS]), 4, [0], err.([BKS]), fmt="s", color="blue")
errorbar(value.([FLAG21]), 5, [0], err.([FLAG21]), fmt="s", color="black")
errorbar(value.([Bali]), 6, [0], err.([Bali]), fmt="s", color="blue")
errorbar(value.([Strass]), 7, [0], err.([Strass]), fmt="s", color="blue")
for i in 1:3:length(sqrt_t0)
        errorbar(value(sqrt_t0[i]), 7+i, 0, err(sqrt_t0[i]), fmt="<", color="deeppink")
end
for i in 2:3:length(sqrt_t0)
        errorbar(value(sqrt_t0[i]), 7+i, 0, err(sqrt_t0[i]), fmt=">", color="deeppink")
end
for i in 3:3:length(sqrt_t0)
        errorbar(value(sqrt_t0[i]), 7+i, 0, err(sqrt_t0[i]), fmt="x", color="deeppink")
end
v = [i for i in 1:1:68]
fill_betweenx(v, value(FLAG21)-err(FLAG21), value(FLAG21)+err(FLAG21), color="lightsteelblue", alpha=0.40)

yticks([collect(1:7); 7 .+ collect(1:3:3*length(sqrt_t0_c))], 
       [L"BMW '12 $[m_{\Omega}]$",
        L"RBC/UKQCD '14 $[m_{\Omega}]$",
        L"QCDSF/UKQCD '15 $[m_O,\;M_V]$",
        L"Bruno et al. '16 $[f_{\pi K}]$",
        "FLAG '21", 
        L"Bali et al. '22 $[m_{\Xi}]$",
        L"Strassberger '23 $[f_{\pi K}]$",

        L"all, FLAG21, $c_{\beta}=1,\;c_{\phi_2}=0.0035$",
        L"all, FLAG21, $c_{\beta}=c_{\phi_2}=0$",
        L"all, FLAG21, $c_{\beta}=10,\;c_{\phi_2}=0$",
        L"all, FLAG21, $c_{\beta}=10,\;c_{\phi_2}=0.0035$",
        L"all, FLAG21, $c_{\beta}=10,\;c_{\phi_2}=10$",

        L"all, FLAG16, $c_{\beta}=1,\;c_{\phi_2}=0.0035$",
        L"all, FLAG16, $c_{\beta}=c_{\phi_2}=0$",
        L"all, FLAG16, $c_{\beta}=10,\;c_{\phi_2}=0$",
        L"all, FLAG16, $c_{\beta}=10,\;c_{\phi_2}=0.0035$",
        L"all, FLAG16, $c_{\beta}=10,\;c_{\phi_2}=10$",

        L"old, FLAG21, $c_{\beta}=1,\;c_{\phi_2}=0.0035$",
        L"old, FLAG21, $c_{\beta}=c_{\phi_2}=0$",
        L"old, FLAG21, $c_{\beta}=10,\;c_{\phi_2}=0$",
        L"old, FLAG21, $c_{\beta}=10,\;c_{\phi_2}=0.0035$",
        L"old, FLAG21, $c_{\beta}=10,\;c_{\phi_2}=10$",

        L"old, FLAG16, $c_{\beta}=1,\;c_{\phi_2}=0.0035$",
        L"old, FLAG16, $c_{\beta}=c_{\phi_2}=0$",
        L"old, FLAG16, $c_{\beta}=10,\;c_{\phi_2}=0$",
        L"old, FLAG16, $c_{\beta}=10,\;c_{\phi_2}=0.0035$",
        L"old, FLAG16, $c_{\beta}=10,\;c_{\phi_2}=10$"
])

xlabel(L"$\sqrt{t_0}\;\rm{[fm]}$")
ax = gca()
ax[:set_xlim]([0.14, 0.154])
tight_layout()

savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0_variations_together.pdf")
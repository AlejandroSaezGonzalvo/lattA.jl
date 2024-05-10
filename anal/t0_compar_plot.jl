using ADerrors, PyPlot
err = ADerrors.err

BKS = uwreal([0.1467,0.00158],1); uwerr(BKS)
Bali = uwreal([0.4098 / sqrt(8), 0.0025 / sqrt(8)],1); uwerr(Bali)
FLAG21 = uwreal([0.14464,0.00087],1); uwerr(FLAG21)
Strass = uwreal([0.1441,0.0013],1); uwerr(Strass)
sqrt_t0_st = uwreal([0.1436,0.0008],""); uwerr(sqrt_t0_st)
sqrt_t0_tm = uwreal([0.1441,0.0011],""); uwerr(sqrt_t0_tm)
sqrt_t0_comb = uwreal([0.1440,0.0007],""); uwerr(sqrt_t0_comb)

fig = figure("",figsize=(5,5))
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 20

#errorbar(value.([sqrt_t0_st_SU2_fpi]), 1, [0], err.([sqrt_t0_st_SU2_fpi]), fmt="s", color="deeppink")
#errorbar(value.([sqrt_t0_st_SU2_fpik]), 3, [0], err.([sqrt_t0_st_SU2_fpik]), fmt="s", color="deeppink")

#errorbar(value.([sqrt_t0_st_old]), 7, [0], err.([sqrt_t0_st_old]), fmt="s", color="green")
errorbar(value.([sqrt_t0_st]), 13, [0], err.([sqrt_t0_st]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_tm]), 17, [0], err.([sqrt_t0_tm]), fmt="s", color="deeppink")
errorbar(value.([sqrt_t0_comb]), 21, [0], err.([sqrt_t0_comb]), fmt="s", color="deeppink")

errorbar(value.([Strass]), 25, [0], err.([Strass]), fmt="s", color="blue")
errorbar(value.([Bali]), 27, [0], err.([Bali]), fmt="s", color="blue")
errorbar(value.([FLAG21]), 29, [0], err.([FLAG21]), fmt="s", color="blue")
errorbar(value.([BKS]), 31, [0], err.([BKS]), fmt="s", color="blue")

xlabel(L"$\sqrt{t_0}\;\rm{[fm]}$")
yticks([13,17,21,25,27,29,31], 
       ["This work, Wilson",
        "This work, Wtm", 
        "This work, combined",  
        "Strassberger '23", 
        "Bali et al.",
        "FLAG '21", 
        "Bruno et al. '16"])
x_plot = [FLAG21 for i in -1:1:32]
v = [i for i in -1:1:32]

fill_betweenx(v, value(FLAG21)-err(FLAG21), value(FLAG21)+err(FLAG21), color="lightsteelblue", alpha=0.20)
ax = gca()
ax[:set_xlim]([0.14, 0.1486])
tight_layout()



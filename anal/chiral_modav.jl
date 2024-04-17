#mods = mods_346
#mods = collect(1:length(sqrt_t0_vec))

#ix = 3
if ix == 3 ## 1: Wtm, 2:Wilson, 3:combined
    mods = mods_c
else
    mods = mods_aux
end

str = ["wtm", "wilson", "combined"]
sqrt_t0_vec = sqrt_t0_ph_vec[ix]
sqrt_t0 = sqrt_t0_ph[ix]
W = W_aux[ix]
pval = pval_aux[ix]
#mods = collect(1:length(pval))

if ix == 3
    for i in 1:length(pval)
        println(string(mods_c[i], " & ", round(pval[i], digits=4), " & ", round(W[i], digits=4), " & ", round(value(sqrt_t0_ph_vec[ix][i]), digits=4), "(", Int(round(err(sqrt_t0_ph_vec[ix][i]) * 1e4, digits=0)), ") \\\\"))
    end
else
    for i in 1:length(pval)
        println(string(mods[i], " & ", round(pval[i], digits=4), " & ", round(W[i], digits=4), " & ", round(value(sqrt_t0_ph_vec[ix][i]), digits=4), "(", Int(round(err(sqrt_t0_ph_vec[ix][i]) * 1e4, digits=0)), ") \\\\"))
    end
end

## plot t0 BMA
fig = figure(str[ix])
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 10

subplot(131) # Create the 2nd axis of a 3x1 arrax of axes
y = collect(1:1:length(mods))
ax = gca() # Get the handle of the current axis
xlabel(L"$\sqrt{t_0}\;\;[fm]$")   
v = [value(sqrt_t0) for i in 1:length(sqrt_t0_vec)]
e = [err(sqrt_t0) for i in 1:length(sqrt_t0_vec)]
fill_betweenx(y, v-e, v+e, color="deepskyblue", alpha=0.75)
errorbar(value.(sqrt_t0_vec), y, 0 .* y, err.(sqrt_t0_vec), fmt="x", color="black")
ax[:set_xlim]([0.139, 0.1477])
ax[:set_ylim]([0, length(mods)+1])
plt.yticks(y, mods)

subplot(132) # Create the 2nd axis of a 3x1 arrax of axes
barh(y, pval, color="deepskyblue")
xlabel("p-value")
ax = gca()
#ax[:set_ylim]([0.0, 0.40])
plt.yticks(y, mods)
#ax.set_yscale("log")
ax[:set_ylim]([0, length(mods)+1])
setp(ax.get_yticklabels(),visible=false)

subplot(133) # Create the 3rd axis of a 3x1 array of axes
ax = gca()
#ax.set_yscale("log") # Set the x axis to a logarithmic scale
barh(y, W, color="deepskyblue")
#ax[:set_ylim]([0.0, 0.1])
xlabel("W")
plt.yticks(y, mods)
#yticks(rotation=90)
ax[:set_ylim]([0, length(mods)+1])
setp(ax.get_yticklabels(),visible=false)

tight_layout()
#subplots_adjust(hspace=0.15) 

#savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/BMA_comb.pdf")


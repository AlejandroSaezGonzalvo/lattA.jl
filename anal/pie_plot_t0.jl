using ADerrors, juobs, PyPlot

fig = figure()

labels = ("Gauge noise", "χ-cont lim. model av.", L"Exp. $f_K$", L"QED $f_K$", L"$|V_{us}|$", L"QED $f_{\pi}$", "Ren. & impr.")
sizes = [53.35, 36, 1.64, 1.64, 6.54, 0.17, 0.66]
explode = [0.05 for i in 1:length(sizes)] # this is the line that separate the slices of the pie

colors= ["red", "orange", "yellow", "green", "blue", "purple", "deeppink"]

pie(sizes,explode=explode,colors=colors)#, autopct="%1.1f%%")
legend(labels, fontsize=10, bbox_to_anchor=(0.5,0.9), loc="upper right")

title(L"$\sqrt{t_0}$ total error squared [Combined]", size=14)
display(gcf())
savefig("t0_combined_pie.pdf")


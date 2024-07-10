using Revise, juobs, BDIO, ADerrors, PyPlot

sqrt_t0_combined_vec = Array{uwreal,1}()
sqrt_t0_vec = Array{uwreal,1}()
sqrt_t0_st_vec = Array{uwreal,1}()

push!(sqrt_t0_combined_vec, sqrt(uwreal([0.1438 ^ 2, 0.0], "t0_guess")))
push!(sqrt_t0_vec, sqrt(uwreal([0.1438 ^ 2, 0.0], "t0_guess")))
push!(sqrt_t0_st_vec, sqrt(uwreal([0.1438 ^ 2, 0.0], "t0_guess")))

fb = BDIO_open("/home/asaez/cls_ens/results/t0_1_newrw.bdio", "r") 
BDIO_seek!(fb) 
aux = read_uwreal(fb)
push!(sqrt_t0_vec, sqrt(aux))
BDIO_seek!(fb, 2)
aux = read_uwreal(fb)
push!(sqrt_t0_st_vec, sqrt(aux))
BDIO_seek!(fb, 2)  
aux = read_uwreal(fb)
push!(sqrt_t0_combined_vec, sqrt(aux))

fb = BDIO_open("/home/asaez/cls_ens/results/t0_2_newrw.bdio", "r") 
BDIO_seek!(fb) 
aux = read_uwreal(fb)
push!(sqrt_t0_vec, sqrt(aux))
BDIO_seek!(fb, 2)
aux = read_uwreal(fb)
push!(sqrt_t0_st_vec, sqrt(aux))
BDIO_seek!(fb, 2)  
aux = read_uwreal(fb)
push!(sqrt_t0_combined_vec, sqrt(aux))

fb = BDIO_open("/home/asaez/cls_ens/results/t0_3_newrw.bdio", "r") 
BDIO_seek!(fb) 
aux = read_uwreal(fb)
push!(sqrt_t0_vec, sqrt(aux))
BDIO_seek!(fb, 2)
aux = read_uwreal(fb)
push!(sqrt_t0_st_vec, sqrt(aux))
BDIO_seek!(fb, 2)  
aux = read_uwreal(fb)
push!(sqrt_t0_combined_vec, sqrt(aux))

fb = BDIO_open("/home/asaez/cls_ens/results/t0_4_newrw.bdio", "r") 
BDIO_seek!(fb) 
aux = read_uwreal(fb)
push!(sqrt_t0_vec, sqrt(aux))
BDIO_seek!(fb, 2)
aux = read_uwreal(fb)
push!(sqrt_t0_st_vec, sqrt(aux))
BDIO_seek!(fb, 2)  
aux = read_uwreal(fb)
push!(sqrt_t0_combined_vec, sqrt(aux))

fb = BDIO_open("/home/asaez/cls_ens/results/t0_5_newrw.bdio", "r") 
BDIO_seek!(fb) 
aux = read_uwreal(fb)
push!(sqrt_t0_vec, sqrt(aux))
BDIO_seek!(fb, 2)
aux = read_uwreal(fb)
push!(sqrt_t0_st_vec, sqrt(aux))
BDIO_seek!(fb, 2)  
aux = read_uwreal(fb)
push!(sqrt_t0_combined_vec, sqrt(aux))

fb = BDIO_open("/home/asaez/cls_ens/results/t0_6_newrw.bdio", "r") 
BDIO_seek!(fb) 
aux = read_uwreal(fb)
push!(sqrt_t0_vec, sqrt(aux))
BDIO_seek!(fb, 2)
aux = read_uwreal(fb)
push!(sqrt_t0_st_vec, sqrt(aux))
BDIO_seek!(fb, 2)  
aux = read_uwreal(fb)
push!(sqrt_t0_combined_vec, sqrt(aux))

uwerr.(sqrt_t0_combined_vec)
uwerr.(sqrt_t0_vec)
uwerr.(sqrt_t0_st_vec)

fig = figure("pyplot_subplot_column10")
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 10
x = collect(0:length(sqrt_t0_combined_vec)-1)
y = sqrt_t0_combined_vec[1:end]
errorbar((x), value.(y), err.(y), fmt="x", label="combined", color="blue")
y = sqrt_t0_vec[1:end]
errorbar((x.+0.2), value.(y), err.(y), fmt="<", label="Wtm", color="darkorange")
y = sqrt_t0_st_vec[1:end]
errorbar((x.+0.1), value.(y), err.(y), fmt=">", label="Wilson", color="green")
xlabel("iteration")
ylabel(L"$\sqrt{t_0}\;\;[fm]$")
xticks(rotation=0)
tight_layout()
legend()  
savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0_iter.pdf")
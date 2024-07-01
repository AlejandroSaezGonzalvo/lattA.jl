using juobs, BDIO, ADerrors

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");

#phi4_ph = uwreal([1.11,0.0],"phi4_ph")
#sqrt_t0_ph = sqrt(phi4_ph / 8 / (MK ^ 2 + 0.5 * Mpi ^ 2) * hc ^ 2)
#t0_ph = sqrt_t0_ph ^ 2

#fb = BDIO_open("/home/asaez/cls_ens/results/t0_2nd.bdio", "r") 
#BDIO_seek!(fb) 
#BDIO_seek!(fb, 2)
#BDIO_seek!(fb, 2)  
#t0_ph = read_uwreal(fb)
#sqrt_t0_ph = sqrt(t0_ph)

sqrt_t0_ph = uwreal([0.1438,0.0],"t0_ph")
t0_ph = sqrt_t0_ph ^ 2

phi4_ph = 8 * sqrt_t0_ph ^ 2 * (MK ^ 2 + 0.5 * Mpi ^ 2) / hc ^ 2
phi2_ph = 8 * sqrt_t0_ph ^ 2 * Mpi ^ 2 / hc ^ 2

fpik_exp = (2/3) * (0.5 * Fpi + FK)
fpik_exp = fpik_exp / hc
uwerr(fpik_exp)

t0_old = t0_ph
LambdaQCD = 340 ## in Mev
Lam = LambdaQCD / hc

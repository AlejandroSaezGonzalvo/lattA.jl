using ADerrors, LsqFit, ForwardDiff, LinearAlgebra, SpecialFunctions

function md_1(x1,x2,p) 
    ## mass derivative for ZA * \sqrt{t0} * f_{\pi K}, \phi_2, \phi_4^{\rm tm}, t0
    return p[1] + p[2] * x1 + p[3] * x2
end

function md_2(x1,x2,p) 
    ## mass derivative for 1 / ZP * m_{12}^{\rm (v)} * \sqrt{t0}
   return p[1] + p[2] * x1 + p[3] * x1 ^ 2 + p[4] * x1 + p[5] * x1 * x2
end

function match_m12_sym(x,p)
    return [p[3] * (1/x[i,1] - 1/p[1]) + p[4] * (x[i,2] - p[2]) for i in 1:length(x[:,1])]
end

function match_phi2_sym(x,p)
    return [p[5]* (1/ x[i,1] - 1/p[1]) ^ 2 / x[i,2] + p[6] * (x[i,2] - p[2]) for i in 1:length(x[:,1])]
end

function match_m12(x,p)
    return [p[4] * (1/x[i,1] - 1/p[1]) + p[5] * (x[i,2] - p[2]) for i in 1:length(x[:,1])]
end

function match_phi2(x,p)
    return [p[6]* (1/ x[i,1] - 1/p[1]) ^ 2 / x[i,2] + p[7] * (x[i,2] - p[2]) for i in 1:length(x[:,1])]
end

function match_phi4(x,p)
    return [p[8]* (1/x[i,1] - 1/p[1]) ^ 2 / x[i,2] + p[9] * (1/x[i,1] - 1/p[1]) ^ 2 / x[i,3] + p[10] * (x[i,2] - p[2]) + p[11] * (x[i,3]-p[3]) for i in 1:length(x[:,1])]
end

function interp_fpik_sym(x,p)
    return [p[1] + p[2] * x[i,2] + p[3] / x[i,1] + p[4] / x[i,1] ^ 2 for i in 1:length(x[:,1])]
end

function interp_fpik_constTr(x,p)
    return [p[1] + p[2] * x[i,2] + p[3] * x[i,3] + p[4] / x[i,1] + p[5] / x[i,1] ^ 2 for i in 1:length(x[:,1])]
end

## chiral-continuum

function model2_ChPT_a1_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[4] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_ChPT_a2_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:L2]
        return [f; g]
    else
        #t0fpik Wilson:
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (-0.111) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[4] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (-0.111) for i in 1:L2]
        return [f; g]
    end
end

function model2_ChPT_a3_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[4] * x[i,1] * x[i,2] for i in 1:L2]
    return [f; g]
end

function model2_ChPT_a4_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[4] * x[i,1] * x[i,2] + p[5] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_ChPT_a1(x,p) 
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,1] for i in 1:L1]
    return f
end

function model2_ChPT_a2(x,p) 
    if x[1,1] == 0
        return [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:length(x[:,1])]
    else
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (-0.111) for i in 1:L1]
        return f
    end
end

function model2_ChPT_a3(x,p) 
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,2] * x[i,1] for i in 1:L1]
    return f
end

function model2_ChPT_a4(x,p) 
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,2] * x[i,1] + p[4] * x[i,1] for i in 1:L1]
    return f
end

function model2_Taylor_a1_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[4] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_Taylor4_a1_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[4] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[5] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_Taylor_a2_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:L2]
        return [f; g]
    else
        #t0fpik Wilson:
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (-0.111) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[4] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (-0.111) for i in 1:L2]
        return [f; g]
    end
end

function model2_Taylor_a4_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[4] * x[i,1] * x[i,2] + p[5] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_Taylor_a1(x,p) 
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i,1] for i in 1:L1]
    return f
end

function model2_Taylor4_a1(x,p) 
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[4] * x[i,1] for i in 1:L1]
    return f
end

function model2_Taylor_a2(x,p) 
    if x[1,1] == 0
        return [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:length(x[:,1])]
    else
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (-0.111) for i in 1:L1]
        return f
    end
end

function model2_Taylor_a4(x,p) 
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i,2] * x[i,1] + p[4] * x[i,1] for i in 1:L1]
    return f
end

function model2_ChPT2_a1(x,p) 
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[5] * x[i,1] for i in 1:div(length(x[:,2]),2)]
    kaon = [p[3] * x[i,2] + p[4] * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[6] * x[i,1] for i in 1:div(length(x[:,2]),2)]
    return [pion;kaon]
end

function model2_ChPT2_a1_combined(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[5] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    kaon = [p[3] * x[i,2] + p[4] * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[6] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[7] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    ka_tm = [p[3] * x[i,2] + p[4] * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[8] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    return [pion; kaon; pi_tm; ka_tm]
end

function model2_ChPT2_fpi_a1(x,p) 
    return [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,1] for i in 1:length(x[:,2])]
end

function model2_ChPT2_fpi_a1_combined(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,1] for i in 1:div(length(x[:,2]),2)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[4] * x[i,1] for i in 1:div(length(x[:,2]),2)]
    return [pion; pi_tm]
end

function model_plot(x,p) 
    return [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,1] for i in 1:length(x[:,1])]
end

function model_plot_SU2_pi(x,p) 
    return [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[5] * x[i,1] for i in 1:length(x[:,2])]
end

function model_plot_SU2_k(x,p) 
    return [p[3] * x[i,2] + p[4] * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[6] * x[i,1] for i in 1:length(x[:,2])]
end
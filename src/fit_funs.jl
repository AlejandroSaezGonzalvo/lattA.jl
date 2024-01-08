using ADerrors, LsqFit, ForwardDiff, LinearAlgebra, SpecialFunctions

function md_1(x1,x2,p) 
    ## mass derivative for ZA * \sqrt{t0} * f_{\pi K}, \phi_2, \phi_4^{\rm tm}, t0
    return p[1] + p[2] * x1 + p[3] * x2
end

function md_2(x1,x2,p) 
    ## mass derivative for 1 / ZP * m_{12}^{\rm (v)} * \sqrt{t0}
   return p[1] + p[2] * x1 + p[3] * x1 ^ 2 + p[4] * x1 + p[5] * x1 * x2
end

function match_sym(x,p)
    ## m12:
    f = [p[3] * (1/x[i,1] - 1/p[1]) + p[4] * (x[i,2] - p[2]) for i in 1:length(x[:,1])]
    ## phi4:
    h = [p[5]* (1/ x[i,1] - 1/p[1]) ^ 2 / x[i,2] + p[6] * (x[i,2] - p[2]) for i in 1:length(x[:,1])]
    return [f;h]
end

function interp_fpik_sym(x,p)
    return [p[1] + p[2] * x[i,2] + p[3] / x[i,1] + p[4] / x[i,1] ^ 2 for i in 1:length(x[:,1])]
end

function interp_fpik(x,p)
    return [p[1] + p[2] * x[i,2] + p[3] * x[i,3] + p[4] / x[i,1] + p[5] / x[i,1] ^ 2 for i in 1:length(x[:,1])]
end
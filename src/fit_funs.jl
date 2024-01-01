using ADerrors, LsqFit, ForwardDiff, LinearAlgebra, SpecialFunctions

function md_1(x1,x2,p) 
    ## mass derivative for ZA * \sqrt{t0} * f_{\pi K}, \phi_2, \phi_4^{\rm tm}, t0
    return p[1] + p[2] * x1 + p[3] * x2
end

function md_2(x1,x2,p) 
    ## mass derivative for 1 / ZP * m_{12}^{\rm (v)} * \sqrt{t0}
   return p[1] + p[2] * x1 + p[3] * x1 ^ 2 + p[4] * x1 + p[5] * x1 * x2
end
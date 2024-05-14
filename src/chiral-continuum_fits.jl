
Gamma_1 = -0.111
Gamma_2 = 0.198
Gamma_3 = 0.247
Gamma_4 = 0.519
Gamma_5 = 0.583

function model2_ChPT_a2_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[4] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_ChPT_a4_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,1] + p[4] * x[i,1] ^ 2 for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[5] * x[i,1] + p[6] * x[i,1] ^ 2 for i in 1:L2]
    return [f; g]
end

function model2_ChPT_aas_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:L2]
        return [f; g]
    else
        #t0fpik Wilson:
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[4] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:L2]
        return [f; g]
    end
end

function model2_ChPT_aas2_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:L2]
        return [f; g]
    else
        #t0fpik Wilson:
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_2) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[4] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_2) for i in 1:L2]
        return [f; g]
    end
end

function model2_ChPT_aas3_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:L2]
        return [f; g]
    else
        #t0fpik Wilson:
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_3) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[4] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_3) for i in 1:L2]
        return [f; g]
    end
end

function model2_ChPT_aas4_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:L2]
        return [f; g]
    else
        #t0fpik Wilson:
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_4) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[4] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_4) for i in 1:L2]
        return [f; g]
    end
end

function model2_ChPT_aas5_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:L2]
        return [f; g]
    else
        #t0fpik Wilson:
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_5) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[4] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_5) for i in 1:L2]
        return [f; g]
    end
end

function model2_ChPT_a2a2phi2_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[4] * x[i,1] * x[i,2] + p[5] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_ChPT_aas5a2phi2_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:L2]
        return [f; g]
    else
        #t0fpik Wilson:
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_5) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[4] * x[i,1] * x[i,2] + p[5] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_5) for i in 1:L2]
        return [f; g]
    end
end


function model2_ChPT_a2phi2a2_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,1] * x[i,2] + p[4] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[5] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_ChPT_a2phi2_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,1] * x[i,2] + p[4] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[5] * x[i,1] * x[i,2] + p[6] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_ChPT_a2(x,p) 
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,1] for i in 1:L1]
    return f
end

function model2_ChPT_a4(x,p) 
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,1] + p[4] * x[i,1] ^ 2 for i in 1:L1]
    return f
end

function model2_ChPT_aas(x,p) 
    if x[1,1] == 0
        return [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:length(x[:,1])]
    else
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:L1]
        return f
    end
end

function model2_ChPT_aas2(x,p) 
    if x[1,1] == 0
        return [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:length(x[:,1])]
    else
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_2) for i in 1:L1]
        return f
    end
end

function model2_ChPT_aas3(x,p) 
    if x[1,1] == 0
        return [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:length(x[:,1])]
    else
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_3) for i in 1:L1]
        return f
    end
end

function model2_ChPT_aas4(x,p) 
    if x[1,1] == 0
        return [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:length(x[:,1])]
    else
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_4) for i in 1:L1]
        return f
    end
end

function model2_ChPT_aas5(x,p) 
    if x[1,1] == 0
        return [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:length(x[:,1])]
    else
        f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_5) for i in 1:L1]
        return f
    end
end

function model2_ChPT_a2phi2(x,p) 
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,2] * x[i,1] + p[4] * x[i,1] for i in 1:L1]
    return f
end

function model2_ChPT_aas5phi2(x,p) 
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,2] * x[i,1] + p[4] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_5) for i in 1:L1]
    return f
end

function model2_Taylor_a2_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[4] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_Taylor4_a2a2phi2_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[4] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[5] * x[i,1] * x[i,2] + p[6] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_Taylor4_a2phi2a2_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[4] * x[i,1] * x[i,2] + p[5] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[6] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_Taylor4_a2phi2_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[4] * x[i,1] * x[i,2] + p[5] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[6] * x[i,1] * x[i,2] + p[7] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_Taylor4_a2_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[4] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[5] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_Taylor4_aas_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2 + p[3] * (x[i,2] - x[i,4]) ^ 4) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2 + p[3] * (x[i,2] - x[i,4]) ^ 4) for i in 1:L2]
        return [f; g]
    else
        #t0fpik Wilson:
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[4] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[5] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:L2]
        return [f; g]
    end
end

function model2_Taylor_aas_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:L2]
        return [f; g]
    else
        #t0fpik Wilson:
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[4] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:L2]
        return [f; g]
    end
end

function model2_Taylor_aas2_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:L2]
        return [f; g]
    else
        #t0fpik Wilson:
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_2) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[4] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_2) for i in 1:L2]
        return [f; g]
    end
end

function model2_Taylor_aas3_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:L2]
        return [f; g]
    else
        #t0fpik Wilson:
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_3) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[4] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_3) for i in 1:L2]
        return [f; g]
    end
end

function model2_Taylor_aas4_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:L2]
        return [f; g]
    else
        #t0fpik Wilson:
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_4) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[4] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_4) for i in 1:L2]
        return [f; g]
    end
end

function model2_Taylor_aas5_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:L2]
        return [f; g]
    else
        #t0fpik Wilson:
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_5) for i in 1:L1]
        #t0fpik Wtm:
        g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[4] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_5) for i in 1:L2]
        return [f; g]
    end
end

function model2_Taylor_a2a2phi2_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[4] * x[i,1] * x[i,2] + p[5] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_Taylor_a2phi2a2_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i,1] * x[i,2] + p[4] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[5] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_Taylor_a2phi2_combined(x,p)
    #t0fpik Wilson:
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i,1] * x[i,2] + p[4] * x[i,1] for i in 1:L1]
    #t0fpik Wtm:
    g = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[5] * x[i,1] * x[i,2] + p[6] * x[i,1] for i in 1:L2]
    return [f; g]
end

function model2_Taylor_a2(x,p) 
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i,1] for i in 1:L1]
    return f
end

function model2_Taylor4_a2(x,p) 
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[4] * x[i,1] for i in 1:L1]
    return f
end

function model2_Taylor4_aas(x,p) 
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[4] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:L1]
    return f
end

function model2_Taylor4_a2phi2(x,p) 
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * (x[i,2] - x[i,4]) ^ 4 + p[4] * x[i,1] * x[i,2] + p[5] * x[i,1] for i in 1:L1]
    return f
end

function model2_Taylor_aas(x,p) 
    if x[1,1] == 0
        return [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:length(x[:,1])]
    else
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:L1]
        return f
    end
end

function model2_Taylor_aas2(x,p) 
    if x[1,1] == 0
        return [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:length(x[:,1])]
    else
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_2) for i in 1:L1]
        return f
    end
end

function model2_Taylor_aas3(x,p) 
    if x[1,1] == 0
        return [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:length(x[:,1])]
    else
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_3) for i in 1:L1]
        return f
    end
end

function model2_Taylor_aas4(x,p) 
    if x[1,1] == 0
        return [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:length(x[:,1])]
    else
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_4) for i in 1:L1]
        return f
    end
end

function model2_Taylor_aas5(x,p) 
    if x[1,1] == 0
        return [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) for i in 1:length(x[:,1])]
    else
        f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_5) for i in 1:L1]
        return f
    end
end

function model2_Taylor_a2phi2(x,p) 
    f = [(p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2) + p[3] * x[i,2] * x[i,1] + p[4] * x[i,1] for i in 1:L1]
    return f
end







function model2_ChPT3_a2(x,p) 
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] for i in 1:div(length(x[:,2]),2)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + p[8] * x[i,1] for i in 1:div(length(x[:,2]),2)]
    return [pion;kaon]
end

function model2_ChPT3_a2_combined(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + p[8] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + p[9] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + p[10] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    return [pion; kaon; pi_tm; ka_tm]
end

function model2_ChPT3_aas(x,p) 
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:div(length(x[:,2]),2)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + p[8] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:div(length(x[:,2]),2)]
    return [pion;kaon]
end

function model2_ChPT3_aas_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] for i in 1:div(length(x[:,2]),4)]
        kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) for i in 1:div(length(x[:,2]),4)]
        #t0fpik Wtm:
        pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] for i in 1:div(length(x[:,2]),4)]
        ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) for i in 1:div(length(x[:,2]),4)]
        return [pion; kaon; pi_tm; ka_tm]
    else
        #t0fpik Wilson:
        pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:div(length(x[:,2]),4)]
        kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + p[8] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:div(length(x[:,2]),4)]
        #t0fpik Wtm:
        pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + p[9] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:div(length(x[:,2]),4)]
        ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + p[10] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:div(length(x[:,2]),4)]
        return [pion; kaon; pi_tm; ka_tm]
    end
end

function model2_ChPT3_a2phi2(x,p) 
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + (p[4] + p[8] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),2)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),2)]
    return [pion;kaon]
end

function model2_ChPT3_a2phi2_sym(x,p) 
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + (p[4] + p[8] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),2)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + (p[4] + p[8] * (x[i,3] - 0.5 * x[i,2])) * x[i,1] for i in 1:div(length(x[:,2]),2)]
    return [pion;kaon]
end

function model2_ChPT3_a2phi2_combined(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + (p[4] + p[8] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + (p[11] + p[12] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + (p[13] + p[14] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    return [pion; kaon; pi_tm; ka_tm]
end

function model2_ChPT3_a2phi2_combined_sym(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + (p[4] + p[8] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + (p[4] + p[8] * (x[i,3] - 0.5 * x[i,2])) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + (p[9] + p[10] * (x[i,3] - 0.5 * x[i,2])) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    return [pion; kaon; pi_tm; ka_tm]
end

function model2_ChPT3_a2a2phi2_combined(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + p[8] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + (p[11] + p[12] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    return [pion; kaon; pi_tm; ka_tm]
end

function model2_ChPT3_a2a2phi2_combined_sym(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + p[8] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + (p[9] + p[10] * (x[i,3] - 0.5 * x[i,2])) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    return [pion; kaon; pi_tm; ka_tm]
end

function model2_ChPT3_a2phi2a2_combined(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + (p[4] + p[8] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + p[11] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + p[12] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    return [pion; kaon; pi_tm; ka_tm]
end

function model2_ChPT3_a2phi2a2_combined_sym(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + (p[4] + p[8] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + (p[4] + p[8] * (x[i,3] - 0.5 * x[i,2])) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + p[9] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + p[10] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    return [pion; kaon; pi_tm; ka_tm]
end








function model2_ChPT2_a2(x,p) 
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] for i in 1:div(length(x[:,2]),2)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[8] * x[i,1] for i in 1:div(length(x[:,2]),2)]
    return [pion;kaon]
end

function model2_ChPT2_a2_combined(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[8] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[9] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[10] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    return [pion; kaon; pi_tm; ka_tm]
end

function model2_ChPT2_aas(x,p) 
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:div(length(x[:,2]),2)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[8] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:div(length(x[:,2]),2)]
    return [pion;kaon]
end

function model2_ChPT2_aas_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] for i in 1:div(length(x[:,2]),4)]
        kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) for i in 1:div(length(x[:,2]),4)]
        #t0fpik Wtm:
        pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] for i in 1:div(length(x[:,2]),4)]
        ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) for i in 1:div(length(x[:,2]),4)]
        return [pion; kaon; pi_tm; ka_tm]
    else
        #t0fpik Wilson:
        pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:div(length(x[:,2]),4)]
        kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[8] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:div(length(x[:,2]),4)]
        #t0fpik Wtm:
        pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[9] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:div(length(x[:,2]),4)]
        ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[10] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:div(length(x[:,2]),4)]
        return [pion; kaon; pi_tm; ka_tm]
    end
end

function model2_ChPT2_a2phi2(x,p) 
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + (p[4] + p[8] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),2)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),2)]
    return [pion;kaon]
end

function model2_ChPT2_a2phi2_sym(x,p) 
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + (p[4] + p[8] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),2)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + (p[4] + p[8] * (x[i,3] - 0.5 * x[i,2])) * x[i,1] for i in 1:div(length(x[:,2]),2)]
    return [pion;kaon]
end

function model2_ChPT2_a2phi2_combined(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + (p[4] + p[8] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + (p[11] + p[12] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + (p[13] + p[14] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    return [pion; kaon; pi_tm; ka_tm]
end

function model2_ChPT2_a2phi2_combined_sym(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + (p[4] + p[8] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + (p[4] + p[8] * (x[i,3] - 0.5 * x[i,2])) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + (p[9] + p[10] * (x[i,3] - 0.5 * x[i,2])) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    return [pion; kaon; pi_tm; ka_tm]
end

function model2_ChPT2_a2a2phi2_combined(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[8] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + (p[11] + p[12] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    return [pion; kaon; pi_tm; ka_tm]
end

function model2_ChPT2_a2a2phi2_combined_sym(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[8] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + (p[9] + p[10] * (x[i,3] - 0.5 * x[i,2])) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    return [pion; kaon; pi_tm; ka_tm]
end

function model2_ChPT2_a2phi2a2_combined(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + (p[4] + p[8] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + (p[9] + p[10] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[11] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[12] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    return [pion; kaon; pi_tm; ka_tm]
end

function model2_ChPT2_a2phi2a2_combined_sym(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + (p[4] + p[8] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    kaon = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + (p[4] + p[8] * (x[i,3] - 0.5 * x[i,2])) * x[i,1] for i in 1:div(length(x[:,2]),4)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[9] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    ka_tm = [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[10] * x[i,1] for i in 1:div(length(x[:,2]),4)]
    return [pion; kaon; pi_tm; ka_tm]
end








function model_ChPT2_a2(x,p) 
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] for i in 1:length(x[:,2])]
    return pion
end

function model_ChPT2_a2_combined(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] for i in 1:div(length(x[:,2]),2)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[5] * x[i,1] for i in 1:div(length(x[:,2]),2)]
    return [pion; pi_tm]
end

function model_ChPT2_aas(x,p) 
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:length(x[:,2])]
    return pion
end

function model_ChPT2_aas_combined(x,p)
    if x[1,1] == 0
        #t0fpik Wilson:
        pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] for i in 1:div(length(x[:,2]),2)]
        #t0fpik Wtm:
        pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] for i in 1:div(length(x[:,2]),2)]
        return [pion; pi_tm]
    else
        #t0fpik Wilson:
        pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:div(length(x[:,2]),2)]
        #t0fpik Wtm:
        pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[5] * x[i,1] * (-1 / (log(sqrt(x[i,1] * 8 * value(t0_old)) * Lam))) ^ (Gamma_1) for i in 1:div(length(x[:,2]),2)]
        return [pion; pi_tm]
    end
end

function model_ChPT2_a2phi2(x,p) 
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + (p[4] + p[5] * x[i,2]) * x[i,1] for i in 1:length(x[:,2])]
    return pion
end

function model_ChPT2_a2phi2_combined(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + (p[4] + p[5] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),2)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + (p[6] + p[7] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),2)]
    return [pion; kaon; pi_tm; ka_tm]
end

function model_ChPT2_a2a2phi2_combined(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] for i in 1:div(length(x[:,2]),2)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + (p[5] + p[6] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),2)]
    return [pion; pi_tm]
end

function model_ChPT2_a2phi2a2_combined(x,p)
    #t0fpik Wilson:
    pion = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + (p[4] + p[5] * x[i,2]) * x[i,1] for i in 1:div(length(x[:,2]),2)]
    #t0fpik Wtm:
    pi_tm = [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[6] * x[i,1] for i in 1:div(length(x[:,2]),2)]
    return [pion; pi_tm]
end









function model_plot(x,p) 
    return [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) + p[3] * x[i,1] for i in 1:length(x[:,1])]
end

function model_plot_SU2_pi(x,p) 
    return [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] for i in 1:length(x[:,2])]
end

function model_plot_SU2_k(x,p) 
    return [p[5] * x[i,2] + p[6] * x[i,3] + p[7] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[8] * x[i,1] for i in 1:length(x[:,2])]
end

function model_plot_SU3_pi(x,p) 
    return [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2)) + p[3] * x[i,3] + p[4] * x[i,1] for i in 1:length(x[:,2])]
end

function model_plot_SU3_k(x,p) 
    return [p[5] * x[i,2] + p[6] * x[i,3] + p[2] / (4 * pi) * (1 - 3/4 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2) - 3/2 * (x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2 * log((x[i,3] - 0.5 * x[i,2]) / p[2] ^ 2) - 3/4 * (4/3 * x[i,3] - x[i,2]) / p[2] ^ 2 * log((4/3 * x[i,3] - x[i,2]) / p[2] ^ 2)) + p[8] * x[i,1] for i in 1:length(x[:,2])]
end

##########################################################
##########################################################
##########################################################

function fun_phi12_X(x,p)
    ## p = 3
    return [p[1] + p[2] * x[i,2] + p[3] * x[i,5] * (x[i,2] * log(x[i,2]) - 1 / 3 * (4 * x[i,3] - 3 * x[i,2]) / 3 * log((4 * x[i,3] - 3 * x[i,2]) / 3)) for i in 1:length(x[:,2])]
end

function fun_phi13_X(x,p)
    ## p = 3
    return [p[1] + p[2] * (x[i,3] - 0.5 * x[i,2]) + 2 / 3 * p[3] * x[i,5] * (4 * x[i,3] - 3 * x[i,2]) / 3 * log((4 * x[i,3] - 3 * x[i,2]) / 3) for i in 1:length(x[:,2])]
end

function fun_y1(x,p)
    ## p = 2
    return [2 * x[i,2] / (2*x[i,3]-x[i,2]) * (1 + p[2]/p[1] * (3/2 * x[i,2] - x[i,3]) - x[i,5] * (x[i,2] * log(x[i,2]) - (4 * x[i,3] - 3 * x[i,2]) / 3 * log((4 * x[i,3] - 3 * x[i,2]) / 3))) for i in 1:length(x[:,2])]
end

function fun_y2(x,p) 
    ## p = 4
    return [3 * p[1] + 2 * p[2] * x[i,3] + p[3] * x[i,5] * (x[i,2] * log(x[i,2]) + (4 * x[i,3] - 3 * x[i,2]) / 3 * log((4 * x[i,3] - 3 * x[i,2]) / 3)) for i in 1:length(x[:,2])]
end

function fun_phi12_T(x,p)
    ## p = 3
    return [p[1] + p[2] * (x[i,2] - x[i,4]) + p[3] * (x[i,2] - x[i,4]) ^ 2 for i in 1:length(x[:,2])]
end

function fun_phi13_T(x,p)
    ## p = 5
    return [p[1] + p[4] * (x[i,2] - x[i,4]) + p[5] * (x[i,2] - x[i,4]) ^ 2 for i in 1:length(x[:,2])]
end

function fun_phi13_T_v2(x,p)
    ## p = 5
    return [p[3] + p[4] * (x[i,2] - x[i,4]) + p[5] * (x[i,2] - x[i,4]) ^ 2 for i in 1:length(x[:,2])]
end

function fun_phi13_T_v3(x,p)
    ## p = 3
    return [p[1] + p[2] * (x[i,2] - x[i,4]) + p[3] * (x[i,2] - x[i,4]) ^ 2 for i in 1:length(x[:,2])]
end


## model 1, p=5

    function phi12_model1(x,p)
        return fun_phi12_X(x,p) .+ [p[4] * x[i,1] for i in 1:length(x[:,1])]
    end

    function phi13_model1(x,p)
        return fun_phi13_X(x,p) .+ [p[5] * x[i,1] for i in 1:length(x[:,1])]
    end

    function x_model1(x,p)
        return [p[i+5] for i in 1:length(x[:,1])]
    end

## model 1.5, p=7

    function phi12_model15(x,p)
        return fun_phi12_X(x,p) .+ [p[4] * x[i,1] + p[5] * x[i,1] * x[i,2] for i in 1:length(x[:,1])]
    end

    function phi13_model15(x,p)
        return fun_phi13_X(x,p) .+ [p[6] * x[i,1] + p[7] * x[i,1] * x[i,2] for i in 1:length(x[:,1])]
    end

    function x_model15(x,p)
        return [p[i+7] for i in 1:length(x[:,1])]
    end

## model 2, p=5

    function y1_model2(x,p)
        return fun_y1(x,p) .+ [p[4] * (2 * x[i,3] - 3 * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
    end

    function phi13_model2(x,p)
        return fun_phi13_X(x,p) .+ [p[5] * x[i,1] for i in 1:length(x[:,1])]
    end

    function x_model2(x,p)
        return [p[i+5] for i in 1:length(x[:,1])]
    end

## model 2.5, p=7

    function y1_model25(x,p)
        return fun_y1(x,p) .+ [(p[4] + p[5] * x[i,2]) * (2 * x[i,3] - 3 * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
    end

    function phi13_model25(x,p)
        return fun_phi13_X(x,p) .+ [(p[6] + p[7] * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
    end

    function x_model25(x,p)
        return [p[i+7] for i in 1:length(x[:,1])]
    end

## model 3, p=5

    function y1_model3(x,p)
        return fun_y1(x,p) .+ [p[4] * (2 * x[i,3] - 3 * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
    end

    function y2_model3(x,p)
        return fun_y2(x,p) .+ [p[5] * x[i,1] for i in 1:length(x[:,1])]
    end

    function x_model3(x,p)
        return [p[i+5] for i in 1:length(x[:,1])]
    end

## model 3.5, p=7

    function y1_model35(x,p)
        return fun_y1(x,p) .+ [(p[4] + p[5] * x[i,2]) * (2 * x[i,3] - 3 * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
    end

    function y2_model35(x,p)
        return fun_y2(x,p) .+ [(p[6] + p[7] * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
    end

    function x_model35(x,p)
        return [p[i+7] for i in 1:length(x[:,1])]
    end

## model 4, p=7

    function y1_model4(x,p)
        return fun_y1(x,p) .+ [p[6] * (2 * x[i,3] - 3 * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
    end

    function phi13_model4(x,p)
        return fun_phi13_T_v2(x,p) .+ [p[7] * x[i,1] for i in 1:length(x[:,1])] 
    end

    function x_model4(x,p)
        return [p[i+7] for i in 1:length(x[:,1])]
    end

## model 4.5, p=9

    function y1_model45(x,p)
        return fun_y1(x,p) .+ [(p[6] + p[7] * x[i,2]) * (2 * x[i,3] - 3 * x[i,2]) * x[i,1] for i in 1:length(x[:,1])]
    end

    function phi13_model45(x,p)
        return fun_phi13_T_v2(x,p) .+ [(p[8] + p[9] * x[i,2]) * x[i,1] for i in 1:length(x[:,1])] 
    end

    function x_model45(x,p)
        return [p[i+9] for i in 1:length(x[:,1])]
    end

## model 5

    function phi12_model5(x,p)
        return fun_phi12_X(x,p) .+ [p[4] * x[i,1] for i in 1:length(x[:,1])]
    end

    function x_model5(x,p)
        return [p[i+5] for i in 1:length(x[:,1])]
    end

## model 6

    function phi13_model6(x,p)
        return fun_phi13_X(x,p) .+ [p[4] * x[i,1] for i in 1:length(x[:,1])]
    end

    function x_model6(x,p)
        return [p[i+5] for i in 1:length(x[:,1])]
    end

## model 7

    function phi12_model7(x,p)
        return fun_phi12_T(x,p) .+ [p[4] * x[i,1] for i in 1:length(x[:,1])] 
    end

    function x_model7(x,p)
        return [p[i+5] for i in 1:length(x[:,1])]
    end

## model 8

    function phi13_model8(x,p)
        return fun_phi13_T_v3(x,p) .+ [p[4] * x[i,1] for i in 1:length(x[:,1])]
    end

    function x_model8(x,p)
        return [p[i+5] for i in 1:length(x[:,1])]
    end








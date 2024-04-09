function matching_sym_plot()
    subplot(221) # Create the 1st axis of a 3x1 array of axes
    ax = gca() # Get the handle of the current axis
    delta = (kappa[2]-kappa[1])/10
    xp = collect(kappa[1]-2*delta:delta:kappa[end]+2*delta) .* [uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,2] = up[2] for i in 1:length(xp[:,1])]
    ylabel(L"$am_{12}^{\rm (v)}$")  
    #xlabel(L"$1/\kappa^{\rm(v)}$")   
    vect = match_m12_sym(xp, up)
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    xp = value.(xp)
    fill_between(1 ./ xp[:,1], v-e, v+e, color="orange", alpha=0.5)
    fill_between(1 ./ xp[:,1], 0, 0, color="gray", alpha=0.5)
    errorbar(1 ./ x[1:3:9,1], value.(m12_sh[1:3:9]), err.(m12_sh[1:3:9]), fmt="x", color="purple", label=L"$a\mu_l^{\rm(v)}=???$")
    errorbar(1 ./ x[2:3:9,1], value.(m12_sh[2:3:9]), err.(m12_sh[2:3:9]), fmt="x", color="green", label=L"$a\mu_l^{\rm(v)}=???$")
    errorbar(1 ./ x[3:3:9,1], value.(m12_sh[3:3:9]), err.(m12_sh[3:3:9]), fmt="x", color="blue", label=L"$a\mu_l^{\rm(v)}=???$")
    kappa_target = 1 / up[1]; uwerr(kappa_target)
    errorbar(value(kappa_target), 0, 0, err(kappa_target), fmt="x", color="black")
    ax[:set_ylim]([value(m12_sh[end])*1.5, value(m12_sh[1])*1.5])
    setp(ax.get_xticklabels(),visible=false)
    #legend()

    subplot(222) 
    ax = gca() 
    delta = (mul[2]-mul[1])/10
    xp = collect(mul[1]-2*delta:delta:mul[end]+2*delta) .* [uwreal([0.0,0.0], "aux") uwreal([1.0,0.0], "aux")]
    [xp[i,1] = up[1] for i in 1:length(xp[:,1])]
    #ylabel(L"$am_{12}^{\rm (v)}$")  
    #xlabel(L"$\mu_l^{\rm(v)}$")   
    vect = match_m12_sym(xp, up)
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    xp = value.(xp)
    fill_between(xp[:,2], v-e, v+e, color="orange", alpha=0.5)
    fill_between(xp[:,2], 0, 0, color="gray", alpha=0.5)
    errorbar(x[1:3,2], value.(m12_sh[1:3]), err.(m12_sh[1:3]), fmt="x", color="purple", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x[4:6,2], value.(m12_sh[4:6]), err.(m12_sh[4:6]), fmt="x", color="green", label=L"$a\kappa^{\rm(v)}=???$")
    errorbar(x[7:9,2], value.(m12_sh[7:9]), err.(m12_sh[7:9]), fmt="x", color="blue", label=L"$\kappa^{\rm(v)}=???$")
    mul_target = up[2]; uwerr(mul_target)
    errorbar(value(mul_target), 0, 0, err(mul_target), fmt="x", color="black")
    ax[:set_ylim]([value(m12_sh[end])*1.5, value(m12_sh[1])*1.5])
    setp(ax.get_xticklabels(),visible=false)
    setp(ax.get_yticklabels(),visible=false)
    
    subplot(223) # Create the 1st axis of a 3x1 array of axes
    ax = gca() # Get the handle of the current axis
    delta = (kappa[2]-kappa[1])/10
    xp = collect(kappa[1]-2*delta:delta:kappa[end]+2*delta) .* [uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,2] = up[2] for i in 1:length(xp[:,1])]
    ylabel(L"$\phi_2^{\rm (v)}$")  
    xlabel(L"$1/\kappa^{\rm(v)}$")   
    vect = match_phi2_sym(xp, up) .+ 2/3 .* [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    phi2_target = 2/3 .* [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(phi2_target)
    v2 = value.(phi2_target)
    e2 = err.(phi2_target)
    xp = value.(xp)
    fill_between(1 ./ xp[:,1], v-e, v+e, color="orange", alpha=0.5)
    fill_between(1 ./ xp[:,1], v2-e2, v2+e2, color="gray", alpha=0.5)
    errorbar(1 ./ x[1:3:9,1], value.(phi2_sh[1:3:9]), err.(phi2_sh[1:3:9]), fmt="x", color="purple", label=L"$a\mu_l^{\rm(v)}=???$")
    errorbar(1 ./ x[2:3:9,1], value.(phi2_sh[2:3:9]), err.(phi2_sh[2:3:9]), fmt="x", color="green", label=L"$a\mu_l^{\rm(v)}=???$")
    errorbar(1 ./ x[3:3:9,1], value.(phi2_sh[3:3:9]), err.(phi2_sh[3:3:9]), fmt="x", color="blue", label=L"$a\mu_l^{\rm(v)}=???$")
    kappa_target = 1 / up[1]; uwerr(kappa_target)
    errorbar(value(kappa_target), value(phi2_target[1]), err(phi2_target[2]), err(kappa_target), fmt="x", color="black")
    ax[:set_ylim]([value(phi2_sh[1])*0.97, value(phi2_sh[end])*1.03])
    xticks(rotation=60)
    #legend()

    subplot(224) 
    ax = gca() 
    delta = (mul[2]-mul[1])/10
    xp = collect(mul[1]-2*delta:delta:mul[end]+2*delta) .* [uwreal([0.0,0.0], "aux") uwreal([1.0,0.0], "aux")]
    [xp[i,1] = up[1] for i in 1:length(xp[:,1])]
    #ylabel(L"$\phi_2^{\rm (v)}$")  
    xlabel(L"$\mu_l^{\rm(v)}$")   
    vect = match_phi2_sym(xp, up) .+ 2/3 .* [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    phi2_target = 2/3 .* [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(phi2_target)
    v2 = value.(phi2_target)
    e2 = err.(phi2_target)
    xp = value.(xp)
    fill_between(xp[:,2], v-e, v+e, color="orange", alpha=0.5)
    fill_between(xp[:,2], v2-e2, v2+e2, color="gray", alpha=0.5)
    errorbar(x[1:3,2], value.(phi2_sh[1:3]), err.(phi2_sh[1:3]), fmt="x", color="purple", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x[4:6,2], value.(phi2_sh[4:6]), err.(phi2_sh[4:6]), fmt="x", color="green", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x[7:9,2], value.(phi2_sh[7:9]), err.(phi2_sh[7:9]), fmt="x", color="blue", label=L"$\kappa^{\rm(v)}=???$")
    mul_target = up[2]; uwerr(mul_target)
    errorbar(value(mul_target), value(phi2_target[1]), err(phi2_target[2]), err(mul_target), fmt="x", color="black")
    ax[:set_ylim]([value(phi2_sh[1])*0.97, value(phi2_sh[end])*1.03])
    setp(ax.get_yticklabels(),visible=false)
    xticks(rotation=60)

    #legend(bbox_to_anchor=(0.5, 0.7))

    tight_layout()
    savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/matching_",id,".pdf"))
    #close("all")
end

function matching_constTr_plot()
    subplot(321) # Create the 1st axis of a 3x1 array of axes
    ax = gca() # Get the handle of the current axis
    delta = (kappa[2]-kappa[1])/10
    xp = collect(kappa[1]-2*delta:delta:kappa[end]+2*delta) .* [uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,2] = up[2] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    ylabel(L"$am_{12}^{\rm (v)}$")  
    #xlabel(L"$1/\kappa^{\rm(v)}$")   
    vect = match_m12(xp, up)
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    xp = value.(xp)
    fill_between(1 ./ xp[:,1], v-e, v+e, color="orange", alpha=0.5)
    fill_between(1 ./ xp[:,1], 0, 0, color="gray", alpha=0.5)
    errorbar(1 ./ x_l[1:3:9,1], value.(m12_sh[1:3:9]), err.(m12_sh[1:3:9]), fmt="x", color="purple", label=L"$a\mu_l^{\rm(v)}=???$") ## mus is fixed
    errorbar(1 ./ x_l[2:3:9,1], value.(m12_sh[2:3:9]), err.(m12_sh[2:3:9]), fmt="x", color="green", label=L"$a\mu_l^{\rm(v)}=???$")
    errorbar(1 ./ x_l[3:3:9,1], value.(m12_sh[3:3:9]), err.(m12_sh[3:3:9]), fmt="x", color="blue", label=L"$a\mu_l^{\rm(v)}=???$")
    kappa_target = 1 / up[1]; uwerr(kappa_target)
    errorbar(value(kappa_target), 0, 0, err(kappa_target), fmt="x", color="black")
    ax[:set_ylim]([value(m12_sh[end])*1.5, value(m12_sh[1])*1.5])
    setp(ax.get_xticklabels(),visible=false)
    #legend()

    subplot(322) 
    ax = gca() 
    delta = (mul[2]-mul[1])/10
    xp = collect(mul[1]-2*delta:delta:mul[end]+2*delta) .* [uwreal([0.0,0.0], "aux") uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,1] = up[1] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    #ylabel(L"$am_{12}^{\rm (v)}$")  
    #xlabel(L"$\mu_l^{\rm(v)}$")   
    vect = match_m12(xp, up)
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    xp = value.(xp)
    fill_between(xp[:,2], v-e, v+e, color="orange", alpha=0.5)
    fill_between(xp[:,2], 0, 0, color="gray", alpha=0.5)
    errorbar(x_l[1:3,2], value.(m12_sh[1:3]), err.(m12_sh[1:3]), fmt="x", color="purple", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x_l[4:6,2], value.(m12_sh[4:6]), err.(m12_sh[4:6]), fmt="x", color="green", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x_l[7:9,2], value.(m12_sh[7:9]), err.(m12_sh[7:9]), fmt="x", color="blue", label=L"$\kappa^{\rm(v)}=???$")
    mul_target = up[2]; uwerr(mul_target)
    errorbar(value(mul_target), 0, 0, err(mul_target), fmt="x", color="black")
    ax[:set_ylim]([value(m12_sh[end])*1.5, value(m12_sh[1])*1.5])
    setp(ax.get_xticklabels(),visible=false)
    setp(ax.get_yticklabels(),visible=false)
    
    subplot(323) # Create the 1st axis of a 3x1 array of axes
    ax = gca() # Get the handle of the current axis
    delta = (kappa[2]-kappa[1])/10
    xp = collect(kappa[1]-2*delta:delta:kappa[end]+2*delta) .* [uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,2] = up[2] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    ylabel(L"$\phi_2^{\rm (v)}$")  
    #xlabel(L"$1/\kappa^{\rm(v)}$")   
    vect =  match_phi2(xp, up) .+ [phi2_w_sh for i in 1:length(xp[:,1])]
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    phi2_target = [phi2_w_sh for i in 1:length(xp[:,1])]
    uwerr.(phi2_target)
    v2 = value.(phi2_target)
    e2 = err.(phi2_target)
    xp = value.(xp)
    fill_between(1 ./ xp[:,1], v-e, v+e, color="orange", alpha=0.5)
    fill_between(1 ./ xp[:,1], v2-e2, v2+e2, color="gray", alpha=0.5)
    errorbar(1 ./ x_l[1:3:9,1], value.(phi2_sh[1:3:9]), err.(phi2_sh[1:3:9]), fmt="x", color="purple", label=L"$a\mu_l^{\rm(v)}=???$") ## mus is fixed
    errorbar(1 ./ x_l[2:3:9,1], value.(phi2_sh[2:3:9]), err.(phi2_sh[2:3:9]), fmt="x", color="green", label=L"$a\mu_l^{\rm(v)}=???$")
    errorbar(1 ./ x_l[3:3:9,1], value.(phi2_sh[3:3:9]), err.(phi2_sh[3:3:9]), fmt="x", color="blue", label=L"$a\mu_l^{\rm(v)}=???$")
    kappa_target = 1 / up[1]; uwerr(kappa_target)
    errorbar(value(kappa_target), value(phi2_target[1]), err(phi2_target[2]), err(kappa_target), fmt="x", color="black")
    ax[:set_ylim]([value(phi2_sh[1])*0.97, value(phi2_sh[end])*1.03])
    setp(ax.get_xticklabels(),visible=false)
    xticks(rotation=60)
    #legend()

    subplot(324) 
    ax = gca() 
    delta = (mul[2]-mul[1])/10
    xp = collect(mul[1]-2*delta:delta:mul[end]+2*delta) .* [uwreal([0.0,0.0], "aux") uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,1] = up[1] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    #ylabel(L"$\phi_2^{\rm (v)}$")  
    #xlabel(L"$\mu_l^{\rm(v)}$")   
    vect =  match_phi2(xp, up) .+ [phi2_w_sh for i in 1:length(xp[:,1])]
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    phi2_target = [phi2_w_sh for i in 1:length(xp[:,1])]
    uwerr.(phi2_target)
    v2 = value.(phi2_target)
    e2 = err.(phi2_target)
    xp = value.(xp)
    fill_between(xp[:,2], v-e, v+e, color="orange", alpha=0.5)
    fill_between(xp[:,2], v2-e2, v2+e2, color="gray", alpha=0.5)
    errorbar(x_l[1:3,2], value.(phi2_sh[1:3]), err.(phi2_sh[1:3]), fmt="x", color="purple", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x_l[4:6,2], value.(phi2_sh[4:6]), err.(phi2_sh[4:6]), fmt="x", color="green", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x_l[7:9,2], value.(phi2_sh[7:9]), err.(phi2_sh[7:9]), fmt="x", color="blue", label=L"$\kappa^{\rm(v)}=???$")
    mul_target = up[2]; uwerr(mul_target)
    errorbar(value(mul_target), value(phi2_target[1]), err(phi2_target[2]), err(mul_target), fmt="x", color="black")
    ax[:set_ylim]([value(phi2_sh[1])*0.97, value(phi2_sh[end])*1.03])
    setp(ax.get_yticklabels(),visible=false)
    setp(ax.get_xticklabels(),visible=false)
    xticks(rotation=60)

    subplot(325) # Create the 1st axis of a 3x1 array of axes
    ax = gca() # Get the handle of the current axis
    delta = (kappa[2]-kappa[1])/10
    xp = collect(kappa[1]-2*delta:delta:kappa[end]+2*delta) .* [uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,2] = up[2] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    ylabel(L"$\phi_4^{\rm (v)}$")  
    xlabel(L"$1/\kappa^{\rm(v)}$")   
    vect = match_phi4(xp, up) .+ [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    phi4_target = [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(phi4_target)
    v2 = value.(phi4_target)
    e2 = err.(phi4_target)
    xp = value.(xp)
    fill_between(1 ./ xp[:,1], v-e, v+e, color="orange", alpha=0.5)
    fill_between(1 ./ xp[:,1], v2-e2, v2+e2, color="gray", alpha=0.5)
    errorbar(1 ./ x_s[1:9:27,1], value.(phi4_sh[1:9:27]), err.(phi4_sh[1:9:27]), fmt="x", color="purple", label=L"$a\mu_l^{\rm(v)}=???$") ## mus is fixed
    errorbar(1 ./ x_s[4:9:27,1], value.(phi4_sh[4:9:27]), err.(phi4_sh[4:9:27]), fmt="x", color="green", label=L"$a\mu_l^{\rm(v)}=???$")
    errorbar(1 ./ x_s[7:9:27,1], value.(phi4_sh[7:9:27]), err.(phi4_sh[7:9:27]), fmt="x", color="blue", label=L"$a\mu_l^{\rm(v)}=???$")
    kappa_target = 1 / up[1]; uwerr(kappa_target)
    errorbar(value(kappa_target), value(phi4_target[1]), err(phi4_target[1]), err(kappa_target), fmt="x", color="black")
    ax[:set_ylim]([value(phi4_sh[1])*0.97, value(phi4_sh[end])*1.03])
    xticks(rotation=60)
    #legend()

    subplot(326) 
    ax = gca() 
    delta = (mul[2]-mul[1])/10
    xp = collect(mul[1]-2*delta:delta:mul[end]+2*delta) .* [uwreal([0.0,0.0], "aux") uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,1] = up[1] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    #ylabel(L"$\phi_2^{\rm (v)}$")  
    xlabel(L"$\mu_l^{\rm(v)}$")   
    vect = match_phi4(xp, up) .+ [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    phi4_target = [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(phi4_target)
    v2 = value.(phi4_target)
    e2 = err.(phi4_target)
    xp = value.(xp)
    fill_between(xp[:,2], v-e, v+e, color="orange", alpha=0.5)
    fill_between(xp[:,2], v2-e2, v2+e2, color="gray", alpha=0.5)
    errorbar(x_s[1:3:9,2], value.(phi4_sh[1:3:9]), err.(phi4_sh[1:3:9]), fmt="x", color="purple", label=L"$\kappa^{\rm(v)}=???$") ## mus is fixed
    errorbar(x_s[10:3:18,2], value.(phi4_sh[10:3:18]), err.(phi4_sh[10:3:18]), fmt="x", color="green", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x_s[19:3:27,2], value.(phi4_sh[19:3:27]), err.(phi4_sh[19:3:27]), fmt="x", color="blue", label=L"$\kappa^{\rm(v)}=???$")
    mul_target = up[2]; uwerr(mul_target)
    errorbar(value(mul_target), value(phi4_target[1]), err(phi4_target[1]), err(mul_target), fmt="x", color="black")
    ax[:set_ylim]([value(phi4_sh[1])*0.97, value(phi4_sh[end])*1.03])
    setp(ax.get_yticklabels(),visible=false)
    xticks(rotation=60)

    #legend(bbox_to_anchor=(0.5, 0.7))

    tight_layout()
    savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/matching_",id,".pdf"))
    #close("all")
end

function matching_D200_J303_plot()
    subplot(221) # Create the 1st axis of a 3x1 array of axes
    ax = gca() # Get the handle of the current axis
    delta = (kappa[2]-kappa[1])/10
    xp = collect(kappa[1]-2*delta:delta:kappa[end]+2*delta) .* [uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,2] = up[2] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    ylabel(L"$am_{12}^{\rm (v)}$")  
    #xlabel(L"$1/\kappa^{\rm(v)}$")   
    vect = match_m12(xp, up)
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    xp = value.(xp)
    fill_between(1 ./ xp[:,1], v-e, v+e, color="orange", alpha=0.5)
    fill_between(1 ./ xp[:,1], 0, 0, color="gray", alpha=0.5)
    errorbar(1 ./ x_l[1:2:6,1], value.(m12_sh[1:2:6]), err.(m12_sh[1:2:6]), fmt="x", color="purple", label=L"$a\mu_l^{\rm(v)}=???$") ## mus is fixed
    errorbar(1 ./ x_l[2:2:6,1], value.(m12_sh[2:2:6]), err.(m12_sh[2:2:6]), fmt="x", color="green", label=L"$a\mu_l^{\rm(v)}=???$")
    kappa_target = 1 / up[1]; uwerr(kappa_target)
    errorbar(value(kappa_target), 0, 0, err(kappa_target), fmt="x", color="black")
    ax[:set_ylim]([value(m12_sh[end])*1.5, value(m12_sh[1])*1.5])
    setp(ax.get_xticklabels(),visible=false)
    #legend()

    subplot(222) 
    ax = gca() 
    delta = (mul[2]-mul[1])/10
    xp = collect(mul[1]-2*delta:delta:mul[end]+2*delta) .* [uwreal([0.0,0.0], "aux") uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,1] = up[1] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    #ylabel(L"$am_{12}^{\rm (v)}$")  
    #xlabel(L"$\mu_l^{\rm(v)}$")   
    vect = match_m12(xp, up)
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    xp = value.(xp)
    fill_between(xp[:,2], v-e, v+e, color="orange", alpha=0.5)
    fill_between(xp[:,2], 0, 0, color="gray", alpha=0.5)
    errorbar(x_l[1:2,2], value.(m12_sh[1:2]), err.(m12_sh[1:2]), fmt="x", color="purple", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x_l[3:4,2], value.(m12_sh[3:4]), err.(m12_sh[3:4]), fmt="x", color="green", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x_l[5:6,2], value.(m12_sh[5:6]), err.(m12_sh[5:6]), fmt="x", color="blue", label=L"$\kappa^{\rm(v)}=???$")
    mul_target = up[2]; uwerr(mul_target)
    errorbar(value(mul_target), 0, 0, err(mul_target), fmt="x", color="black")
    ax[:set_ylim]([value(m12_sh[end])*1.5, value(m12_sh[1])*1.5])
    setp(ax.get_xticklabels(),visible=false)
    setp(ax.get_yticklabels(),visible=false)
    
    subplot(223) # Create the 1st axis of a 3x1 array of axes
    ax = gca() # Get the handle of the current axis
    delta = (kappa[2]-kappa[1])/10
    xp = collect(kappa[1]-2*delta:delta:kappa[end]+2*delta) .* [uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,2] = up[2] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    ylabel(L"$\phi_2^{\rm (v)}$")  
    xlabel(L"$1/\kappa^{\rm(v)}$")   
    vect =  match_phi2(xp, up) .+ [phi2_w_sh for i in 1:length(xp[:,1])]
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    phi2_target = [phi2_w_sh for i in 1:length(xp[:,1])]
    uwerr.(phi2_target)
    v2 = value.(phi2_target)
    e2 = err.(phi2_target)
    xp = value.(xp)
    fill_between(1 ./ xp[:,1], v-e, v+e, color="orange", alpha=0.5)
    fill_between(1 ./ xp[:,1], v2-e2, v2+e2, color="gray", alpha=0.5)
    errorbar(1 ./ x_l[1:2:6,1], value.(phi2_sh[1:2:6]), err.(phi2_sh[1:2:6]), fmt="x", color="purple", label=L"$a\mu_l^{\rm(v)}=???$") ## mus is fixed
    errorbar(1 ./ x_l[2:2:6,1], value.(phi2_sh[2:2:6]), err.(phi2_sh[2:2:6]), fmt="x", color="green", label=L"$a\mu_l^{\rm(v)}=???$")
    kappa_target = 1 / up[1]; uwerr(kappa_target)
    errorbar(value(kappa_target), value(phi2_target[1]), err(phi2_target[2]), err(kappa_target), fmt="x", color="black")
    ax[:set_ylim]([value(phi2_sh[1])*0.97, value(phi2_sh[end])*1.03])
    xticks(rotation=60)
    #legend()

    subplot(224) 
    ax = gca() 
    delta = (mul[2]-mul[1])/10
    xp = collect(mul[1]-2*delta:delta:mul[end]+2*delta) .* [uwreal([0.0,0.0], "aux") uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,1] = up[1] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    #ylabel(L"$\phi_2^{\rm (v)}$")  
    xlabel(L"$\mu_l^{\rm(v)}$")   
    vect =  match_phi2(xp, up) .+ [phi2_w_sh for i in 1:length(xp[:,1])]
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    phi2_target = [phi2_w_sh for i in 1:length(xp[:,1])]
    uwerr.(phi2_target)
    v2 = value.(phi2_target)
    e2 = err.(phi2_target)
    xp = value.(xp)
    fill_between(xp[:,2], v-e, v+e, color="orange", alpha=0.5)
    fill_between(xp[:,2], v2-e2, v2+e2, color="gray", alpha=0.5)
    errorbar(x_l[1:2,2], value.(phi2_sh[1:2]), err.(phi2_sh[1:2]), fmt="x", color="purple", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x_l[3:4,2], value.(phi2_sh[3:4]), err.(phi2_sh[3:4]), fmt="x", color="green", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x_l[5:6,2], value.(phi2_sh[5:6]), err.(phi2_sh[5:6]), fmt="x", color="blue", label=L"$\kappa^{\rm(v)}=???$")
    mul_target = up[2]; uwerr(mul_target)
    errorbar(value(mul_target), value(phi2_target[1]), err(phi2_target[2]), err(mul_target), fmt="x", color="black")
    ax[:set_ylim]([value(phi2_sh[1])*0.97, value(phi2_sh[end])*1.03])
    setp(ax.get_yticklabels(),visible=false)
    xticks(rotation=60)

    #legend(bbox_to_anchor=(0.5, 0.7))

    tight_layout()
    savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/matching_",id,".pdf"))
    #close("all")
end

function matching_newens_sym_plot()
    subplot(221) # Create the 1st axis of a 3x1 array of axes
    ax = gca() # Get the handle of the current axis
    delta = (kappa[2]-kappa[1])/10
    xp = collect(kappa[1]-2*delta:delta:kappa[end]+2*delta) .* [uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,2] = up[2] for i in 1:length(xp[:,1])]
    ylabel(L"$am_{12}^{\rm (v)}$")  
    #xlabel(L"$1/\kappa^{\rm(v)}$")   
    vect = match_m12_sym(xp, up)
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    xp = value.(xp)
    fill_between(1 ./ xp[:,1], v-e, v+e, color="orange", alpha=0.5)
    fill_between(1 ./ xp[:,1], 0, 0, color="gray", alpha=0.5)
    errorbar(1 ./ x[1:2:4,1], value.(m12_sh[1:2:4]), err.(m12_sh[1:2:4]), fmt="x", color="purple", label=L"$a\mu_l^{\rm(v)}=???$")
    errorbar(1 ./ x[2:2:4,1], value.(m12_sh[2:2:4]), err.(m12_sh[2:2:4]), fmt="x", color="green", label=L"$a\mu_l^{\rm(v)}=???$")
    kappa_target = 1 / up[1]; uwerr(kappa_target)
    errorbar(value(kappa_target), 0, 0, err(kappa_target), fmt="x", color="black")
    ax[:set_ylim]([value(m12_sh[end])*1.5, value(m12_sh[1])*1.5])
    setp(ax.get_xticklabels(),visible=false)
    #legend()

    subplot(222) 
    ax = gca() 
    delta = (mul[2]-mul[1])/10
    xp = collect(mul[1]-2*delta:delta:mul[end]+2*delta) .* [uwreal([0.0,0.0], "aux") uwreal([1.0,0.0], "aux")]
    [xp[i,1] = up[1] for i in 1:length(xp[:,1])]
    #ylabel(L"$am_{12}^{\rm (v)}$")  
    #xlabel(L"$\mu_l^{\rm(v)}$")   
    vect = match_m12_sym(xp, up)
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    xp = value.(xp)
    fill_between(xp[:,2], v-e, v+e, color="orange", alpha=0.5)
    fill_between(xp[:,2], 0, 0, color="gray", alpha=0.5)
    errorbar(x[1:2,2], value.(m12_sh[1:2]), err.(m12_sh[1:2]), fmt="x", color="purple", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x[3:4,2], value.(m12_sh[3:4]), err.(m12_sh[3:4]), fmt="x", color="green", label=L"$a\kappa^{\rm(v)}=???$")
    mul_target = up[2]; uwerr(mul_target)
    errorbar(value(mul_target), 0, 0, err(mul_target), fmt="x", color="black")
    ax[:set_ylim]([value(m12_sh[end])*1.5, value(m12_sh[1])*1.5])
    setp(ax.get_xticklabels(),visible=false)
    setp(ax.get_yticklabels(),visible=false)
    
    subplot(223) # Create the 1st axis of a 3x1 array of axes
    ax = gca() # Get the handle of the current axis
    delta = (kappa[2]-kappa[1])/10
    xp = collect(kappa[1]-2*delta:delta:kappa[end]+2*delta) .* [uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,2] = up[2] for i in 1:length(xp[:,1])]
    ylabel(L"$\phi_2^{\rm (v)}$")  
    xlabel(L"$1/\kappa^{\rm(v)}$")   
    vect = match_phi2_sym(xp, up) .+ 2/3 .* [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    phi2_target = 2/3 .* [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(phi2_target)
    v2 = value.(phi2_target)
    e2 = err.(phi2_target)
    xp = value.(xp)
    fill_between(1 ./ xp[:,1], v-e, v+e, color="orange", alpha=0.5)
    fill_between(1 ./ xp[:,1], v2-e2, v2+e2, color="gray", alpha=0.5)
    errorbar(1 ./ x[1:2:4,1], value.(phi2_sh[1:2:4]), err.(phi2_sh[1:2:4]), fmt="x", color="purple", label=L"$a\mu_l^{\rm(v)}=???$")
    errorbar(1 ./ x[2:2:4,1], value.(phi2_sh[2:2:4]), err.(phi2_sh[2:2:4]), fmt="x", color="green", label=L"$a\mu_l^{\rm(v)}=???$")
    kappa_target = 1 / up[1]; uwerr(kappa_target)
    errorbar(value(kappa_target), value(phi2_target[1]), err(phi2_target[2]), err(kappa_target), fmt="x", color="black")
    ax[:set_ylim]([value(phi2_sh[1])*0.97, value(phi2_sh[end])*1.03])
    xticks(rotation=60)
    #legend()

    subplot(224) 
    ax = gca() 
    delta = (mul[2]-mul[1])/10
    xp = collect(mul[1]-2*delta:delta:mul[end]+2*delta) .* [uwreal([0.0,0.0], "aux") uwreal([1.0,0.0], "aux")]
    [xp[i,1] = up[1] for i in 1:length(xp[:,1])]
    #ylabel(L"$\phi_2^{\rm (v)}$")  
    xlabel(L"$\mu_l^{\rm(v)}$")   
    vect = match_phi2_sym(xp, up) .+ 2/3 .* [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    phi2_target = 2/3 .* [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(phi2_target)
    v2 = value.(phi2_target)
    e2 = err.(phi2_target)
    xp = value.(xp)
    fill_between(xp[:,2], v-e, v+e, color="orange", alpha=0.5)
    fill_between(xp[:,2], v2-e2, v2+e2, color="gray", alpha=0.5)
    errorbar(x[1:2,2], value.(phi2_sh[1:2]), err.(phi2_sh[1:2]), fmt="x", color="purple", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x[3:4,2], value.(phi2_sh[3:4]), err.(phi2_sh[3:4]), fmt="x", color="green", label=L"$\kappa^{\rm(v)}=???$")
    mul_target = up[2]; uwerr(mul_target)
    errorbar(value(mul_target), value(phi2_target[1]), err(phi2_target[2]), err(mul_target), fmt="x", color="black")
    ax[:set_ylim]([value(phi2_sh[1])*0.97, value(phi2_sh[end])*1.03])
    setp(ax.get_yticklabels(),visible=false)
    xticks(rotation=60)

    #legend(bbox_to_anchor=(0.5, 0.7))

    tight_layout()
    savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/matching_",id,".pdf"))
    #close("all")
end

function matching_newens_constTr_plot()
    subplot(321) # Create the 1st axis of a 3x1 array of axes
    ax = gca() # Get the handle of the current axis
    delta = (kappa[2]-kappa[1])/10
    xp = collect(kappa[1]-2*delta:delta:kappa[end]+2*delta) .* [uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,2] = up[2] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    ylabel(L"$am_{12}^{\rm (v)}$")  
    #xlabel(L"$1/\kappa^{\rm(v)}$")   
    vect = match_m12(xp, up)
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    xp = value.(xp)
    fill_between(1 ./ xp[:,1], v-e, v+e, color="orange", alpha=0.5)
    fill_between(1 ./ xp[:,1], 0, 0, color="gray", alpha=0.5)
    errorbar(1 ./ x_l[1:2:4,1], value.(m12_sh[1:2:4]), err.(m12_sh[1:2:4]), fmt="x", color="purple", label=L"$a\mu_l^{\rm(v)}=???$") ## mus is fixed
    errorbar(1 ./ x_l[2:2:4,1], value.(m12_sh[2:2:4]), err.(m12_sh[2:2:4]), fmt="x", color="green", label=L"$a\mu_l^{\rm(v)}=???$")
    kappa_target = 1 / up[1]; uwerr(kappa_target)
    errorbar(value(kappa_target), 0, 0, err(kappa_target), fmt="x", color="black")
    ax[:set_ylim]([value(m12_sh[end])*1.5, value(m12_sh[1])*1.5])
    setp(ax.get_xticklabels(),visible=false)
    #legend()

    subplot(322) 
    ax = gca() 
    delta = (mul[2]-mul[1])/10
    xp = collect(mul[1]-2*delta:delta:mul[end]+2*delta) .* [uwreal([0.0,0.0], "aux") uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,1] = up[1] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    #ylabel(L"$am_{12}^{\rm (v)}$")  
    #xlabel(L"$\mu_l^{\rm(v)}$")   
    vect = match_m12(xp, up)
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    xp = value.(xp)
    fill_between(xp[:,2], v-e, v+e, color="orange", alpha=0.5)
    fill_between(xp[:,2], 0, 0, color="gray", alpha=0.5)
    errorbar(x_l[1:2,2], value.(m12_sh[1:2]), err.(m12_sh[1:2]), fmt="x", color="purple", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x_l[3:4,2], value.(m12_sh[3:4]), err.(m12_sh[3:4]), fmt="x", color="green", label=L"$\kappa^{\rm(v)}=???$")
    mul_target = up[2]; uwerr(mul_target)
    errorbar(value(mul_target), 0, 0, err(mul_target), fmt="x", color="black")
    ax[:set_ylim]([value(m12_sh[end])*1.5, value(m12_sh[1])*1.5])
    setp(ax.get_xticklabels(),visible=false)
    setp(ax.get_yticklabels(),visible=false)
    
    subplot(323) # Create the 1st axis of a 3x1 array of axes
    ax = gca() # Get the handle of the current axis
    delta = (kappa[2]-kappa[1])/10
    xp = collect(kappa[1]-2*delta:delta:kappa[end]+2*delta) .* [uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,2] = up[2] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    ylabel(L"$\phi_2^{\rm (v)}$")  
    #xlabel(L"$1/\kappa^{\rm(v)}$")   
    vect =  match_phi2(xp, up) .+ [phi2_w_sh for i in 1:length(xp[:,1])]
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    phi2_target = [phi2_w_sh for i in 1:length(xp[:,1])]
    uwerr.(phi2_target)
    v2 = value.(phi2_target)
    e2 = err.(phi2_target)
    xp = value.(xp)
    fill_between(1 ./ xp[:,1], v-e, v+e, color="orange", alpha=0.5)
    fill_between(1 ./ xp[:,1], v2-e2, v2+e2, color="gray", alpha=0.5)
    errorbar(1 ./ x_l[1:2:4,1], value.(phi2_sh[1:2:4]), err.(phi2_sh[1:2:4]), fmt="x", color="purple", label=L"$a\mu_l^{\rm(v)}=???$") ## mus is fixed
    errorbar(1 ./ x_l[2:2:4,1], value.(phi2_sh[2:2:4]), err.(phi2_sh[2:2:4]), fmt="x", color="green", label=L"$a\mu_l^{\rm(v)}=???$")
    kappa_target = 1 / up[1]; uwerr(kappa_target)
    errorbar(value(kappa_target), value(phi2_target[1]), err(phi2_target[2]), err(kappa_target), fmt="x", color="black")
    ax[:set_ylim]([value(phi2_sh[1])*0.97, value(phi2_sh[end])*1.03])
    setp(ax.get_xticklabels(),visible=false)
    xticks(rotation=60)
    #legend()

    subplot(324) 
    ax = gca() 
    delta = (mul[2]-mul[1])/10
    xp = collect(mul[1]-2*delta:delta:mul[end]+2*delta) .* [uwreal([0.0,0.0], "aux") uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,1] = up[1] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    #ylabel(L"$\phi_2^{\rm (v)}$")  
    #xlabel(L"$\mu_l^{\rm(v)}$")   
    vect =  match_phi2(xp, up) .+ [phi2_w_sh for i in 1:length(xp[:,1])]
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    phi2_target = [phi2_w_sh for i in 1:length(xp[:,1])]
    uwerr.(phi2_target)
    v2 = value.(phi2_target)
    e2 = err.(phi2_target)
    xp = value.(xp)
    fill_between(xp[:,2], v-e, v+e, color="orange", alpha=0.5)
    fill_between(xp[:,2], v2-e2, v2+e2, color="gray", alpha=0.5)
    errorbar(x_l[1:2,2], value.(phi2_sh[1:2]), err.(phi2_sh[1:2]), fmt="x", color="purple", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x_l[3:4,2], value.(phi2_sh[3:4]), err.(phi2_sh[3:4]), fmt="x", color="green", label=L"$\kappa^{\rm(v)}=???$")
    mul_target = up[2]; uwerr(mul_target)
    errorbar(value(mul_target), value(phi2_target[1]), err(phi2_target[2]), err(mul_target), fmt="x", color="black")
    ax[:set_ylim]([value(phi2_sh[1])*0.97, value(phi2_sh[end])*1.03])
    setp(ax.get_xticklabels(),visible=false)
    setp(ax.get_yticklabels(),visible=false)
    xticks(rotation=60)

    subplot(325) # Create the 1st axis of a 3x1 array of axes
    ax = gca() # Get the handle of the current axis
    delta = (kappa[2]-kappa[1])/10
    xp = collect(kappa[1]-2*delta:delta:kappa[end]+2*delta) .* [uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,2] = up[2] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    ylabel(L"$\phi_4^{\rm (v)}$")  
    xlabel(L"$1/\kappa^{\rm(v)}$")   
    vect = match_phi4(xp, up) .+ [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    phi4_target = [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(phi4_target)
    v2 = value.(phi4_target)
    e2 = err.(phi4_target)
    xp = value.(xp)
    fill_between(1 ./ xp[:,1], v-e, v+e, color="orange", alpha=0.5)
    fill_between(1 ./ xp[:,1], v2-e2, v2+e2, color="gray", alpha=0.5)
    errorbar(1 ./ x_s[1:4:8,1], value.(phi4_sh[1:4:8]), err.(phi4_sh[1:4:8]), fmt="x", color="purple", label=L"$a\mu_l^{\rm(v)}=???$") ## mus is fixed
    errorbar(1 ./ x_s[3:4:8,1], value.(phi4_sh[3:4:8]), err.(phi4_sh[3:4:8]), fmt="x", color="green", label=L"$a\mu_l^{\rm(v)}=???$")
    kappa_target = 1 / up[1]; uwerr(kappa_target)
    errorbar(value(kappa_target), value(phi4_target[1]), err(phi4_target[1]), err(kappa_target), fmt="x", color="black")
    ax[:set_ylim]([value(phi4_sh[1])*0.97, value(phi4_sh[end])*1.03])
    xticks(rotation=60)
    #legend()

    subplot(326) 
    ax = gca() 
    delta = (mul[2]-mul[1])/10
    xp = collect(mul[1]-2*delta:delta:mul[end]+2*delta) .* [uwreal([0.0,0.0], "aux") uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,1] = up[1] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    #ylabel(L"$\phi_2^{\rm (v)}$")  
    xlabel(L"$\mu_l^{\rm(v)}$")   
    vect = match_phi4(xp, up) .+ [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    phi4_target = [phi4_ph for i in 1:length(xp[:,1])]
    uwerr.(phi4_target)
    v2 = value.(phi4_target)
    e2 = err.(phi4_target)
    xp = value.(xp)
    fill_between(xp[:,2], v-e, v+e, color="orange", alpha=0.5)
    fill_between(xp[:,2], v2-e2, v2+e2, color="gray", alpha=0.5)
    errorbar(x_s[1:2:4,2], value.(phi4_sh[1:2:4]), err.(phi4_sh[1:2:4]), fmt="x", color="purple", label=L"$\kappa^{\rm(v)}=???$") ## mus is fixed
    errorbar(x_s[5:2:8,2], value.(phi4_sh[5:2:8]), err.(phi4_sh[5:2:8]), fmt="x", color="green", label=L"$\kappa^{\rm(v)}=???$")
    mul_target = up[2]; uwerr(mul_target)
    errorbar(value(mul_target), value(phi4_target[1]), err(phi4_target[1]), err(mul_target), fmt="x", color="black")
    ax[:set_ylim]([value(phi4_sh[1])*0.97, value(phi4_sh[end])*1.03])
    setp(ax.get_yticklabels(),visible=false)
    xticks(rotation=60)

    #legend(bbox_to_anchor=(0.5, 0.7))

    tight_layout()
    savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/matching_",id,".pdf"))
    #close("all")
end

function interp_fpik_sym_plot()
    subplot(121) # Create the 1st axis of a 3x1 array of axes
    ax = gca() # Get the handle of the current axis
    delta = (kappa[2]-kappa[1])/10
    xp = collect(kappa[1]-2*delta:delta:kappa[end]+2*delta) .* [uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,2] = up[2] for i in 1:length(xp[:,1])]
    ylabel(L"$\sqrt{t_0}f_{\pi K}^{\rm (v)}$")  
    xlabel(L"$1/\kappa^{\rm(v)}$")   
    vect = interp_fpik_sym(xp, up_fpik) 
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    xp = value.(xp)
    fill_between(1 ./ xp[:,1], v-e, v+e, color="orange", alpha=0.5)
    errorbar(1 ./ x[1:3:9,1], value.(fpik_sh[1:3:9]), err.(fpik_sh[1:3:9]), fmt="x", color="purple", label=L"$a\mu_l^{\rm(v)}=???$")
    errorbar(1 ./ x[2:3:9,1], value.(fpik_sh[2:3:9]), err.(fpik_sh[2:3:9]), fmt="x", color="green", label=L"$a\mu_l^{\rm(v)}=???$")
    errorbar(1 ./ x[3:3:9,1], value.(fpik_sh[3:3:9]), err.(fpik_sh[3:3:9]), fmt="x", color="blue", label=L"$a\mu_l^{\rm(v)}=???$")
    kappa_target = 1 / up[1]; uwerr(kappa_target)
    errorbar(value(kappa_target), value(fpik_matched), err(fpik_matched), err(kappa_target), fmt="x", color="black")
    ax[:set_ylim]([value(fpik_sh[1])*0.97, value(fpik_sh[end])*1.03])
    xticks(rotation=60)
    #legend()

    subplot(122) 
    ax = gca() 
    delta = (mul[2]-mul[1])/10
    xp = collect(mul[1]-2*delta:delta:mul[end]+2*delta) .* [uwreal([0.0,0.0], "aux") uwreal([1.0,0.0], "aux")]
    [xp[i,1] = up[1] for i in 1:length(xp[:,1])]
    #ylabel(L"$\phi_2^{\rm (v)}$")  
    xlabel(L"$\mu_l^{\rm(v)}$")   
    vect = interp_fpik_sym(xp, up_fpik) 
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    xp = value.(xp)
    fill_between(xp[:,2], v-e, v+e, color="orange", alpha=0.5)
    errorbar(x[1:3,2], value.(fpik_sh[1:3]), err.(fpik_sh[1:3]), fmt="x", color="purple", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x[4:6,2], value.(fpik_sh[4:6]), err.(fpik_sh[4:6]), fmt="x", color="green", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x[7:9,2], value.(fpik_sh[7:9]), err.(fpik_sh[7:9]), fmt="x", color="blue", label=L"$\kappa^{\rm(v)}=???$")
    mul_target = up[2]; uwerr(mul_target)
    errorbar(value(mul_target), value(fpik_matched), err(fpik_matched), err(mul_target), fmt="x", color="black")
    ax[:set_ylim]([value(fpik_sh[1])*0.97, value(fpik_sh[end])*1.03])
    setp(ax.get_yticklabels(),visible=false)
    xticks(rotation=60)

    #legend(bbox_to_anchor=(0.5, 0.7))

    tight_layout()
    savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/interp_fpik_",id,".pdf"))
    #close("all")
end

function interp_fpik_constTr_plot()
    subplot(121) # Create the 1st axis of a 3x1 array of axes
    ax = gca() # Get the handle of the current axis
    delta = (kappa[2]-kappa[1])/10
    xp = collect(kappa[1]-2*delta:delta:kappa[end]+2*delta) .* [uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,2] = up[2] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    ylabel(L"$\sqrt{t_0}f_{\pi K}^{\rm (v)}$")  
    xlabel(L"$1/\kappa^{\rm(v)}$")   
    vect = interp_fpik_constTr(xp, up_fpik) 
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    xp = value.(xp)
    fill_between(1 ./ xp[:,1], v-e, v+e, color="orange", alpha=0.5)
    errorbar(1 ./ x_s[1:9:27,1], value.(fpik_sh[1:9:27]), err.(fpik_sh[1:9:27]), fmt="x", color="purple", label=L"$a\mu_l^{\rm(v)}=???$") ## mus is fixed
    errorbar(1 ./ x_s[4:9:27,1], value.(fpik_sh[4:9:27]), err.(fpik_sh[4:9:27]), fmt="x", color="green", label=L"$a\mu_l^{\rm(v)}=???$")
    errorbar(1 ./ x_s[7:9:27,1], value.(fpik_sh[7:9:27]), err.(fpik_sh[7:9:27]), fmt="x", color="blue", label=L"$a\mu_l^{\rm(v)}=???$")
    kappa_target = 1 / up[1]; uwerr(kappa_target)
    errorbar(value(kappa_target), value(fpik_matched), err(fpik_matched), err(kappa_target), fmt="x", color="black")
    ax[:set_ylim]([value(fpik_sh[1])*0.97, value(fpik_sh[end])*1.03])
    xticks(rotation=60)
    #legend()

    subplot(122) 
    ax = gca() 
    delta = (mul[2]-mul[1])/10
    xp = collect(mul[1]-2*delta:delta:mul[end]+2*delta) .* [uwreal([0.0,0.0], "aux") uwreal([1.0,0.0], "aux") uwreal([0.0,0.0], "aux")]
    [xp[i,1] = up[1] for i in 1:length(xp[:,1])]
    [xp[i,3] = up[3] for i in 1:length(xp[:,1])]
    #ylabel(L"$\phi_2^{\rm (v)}$")  
    xlabel(L"$\mu_l^{\rm(v)}$")   
    vect = interp_fpik_constTr(xp, up_fpik) 
    uwerr.(vect)
    v = value.(vect)
    e = err.(vect)
    xp = value.(xp)
    fill_between(xp[:,2], v-e, v+e, color="orange", alpha=0.5)
    errorbar(x_s[1:3:9,2], value.(fpik_sh[1:3:9]), err.(fpik_sh[1:3:9]), fmt="x", color="purple", label=L"$\kappa^{\rm(v)}=???$") ## mus is fixed
    errorbar(x_s[10:3:18,2], value.(fpik_sh[10:3:18]), err.(fpik_sh[10:3:18]), fmt="x", color="green", label=L"$\kappa^{\rm(v)}=???$")
    errorbar(x_s[19:3:27,2], value.(fpik_sh[19:3:27]), err.(fpik_sh[19:3:27]), fmt="x", color="blue", label=L"$\kappa^{\rm(v)}=???$")
    mul_target = up[2]; uwerr(mul_target)
    errorbar(value(mul_target), value(fpik_matched), err(fpik_matched), err(mul_target), fmt="x", color="black")
    ax[:set_ylim]([value(fpik_sh[1])*0.97, value(fpik_sh[end])*1.03])
    setp(ax.get_yticklabels(),visible=false)
    xticks(rotation=60)

    #legend(bbox_to_anchor=(0.5, 0.7))

    tight_layout()
    savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/interp_fpik_",id,".pdf"))
    #close("all")
end

function t0_model_av(ix::Int64) ## ix -> 1: Wtm, 2:Wilson, 3:combined
    #mods = mods_346

    str = ["wtm", "wilson", "combined"]
    sqrt_t0_vec = sqrt_t0_ph_vec[ix]
    sqrt_t0 = sqrt_t0_ph[ix]
    W = W_aux[ix]
    pval = pval_aux[ix]

    fig = figure(str[ix],figsize=(5,5))
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 12

    subplot(311) 
    #mods = collect(1:length(sqrt_t0_vec))
    x = collect(1:1:length(mods))
    ax = gca() 
    ylabel(L"$\sqrt{t_0}\;\;[fm]$")   
    v = [value(sqrt_t0) for i in 1:length(sqrt_t0_vec)]
    e = [err(sqrt_t0) for i in 1:length(sqrt_t0_vec)]
    fill_between(x, v-e, v+e, color="deepskyblue", alpha=0.75)
    errorbar(x, value.(sqrt_t0_vec), err.(sqrt_t0_vec), fmt="x", color="black")
    ax[:set_ylim]([0.139, 0.1477])
    ax[:set_xlim]([0, length(mods)+1])
    plt.xticks(x, mods)
    setp(ax.get_xticklabels(),visible=false)

    subplot(312) 
    bar(x, pval, color="deepskyblue")
    ylabel("p-value")
    ax = gca()
    #ax[:set_ylim]([0.0, 0.40])
    plt.xticks(x, mods)
    #ax.set_yscale("log")
    ax[:set_xlim]([0, length(mods)+1])
    setp(ax.get_xticklabels(),visible=false)

    subplot(313) 
    ax = gca()
    bar(x, W, color="deepskyblue")
    #ax[:set_ylim]([0.0, 0.1])
    ylabel("W")
    plt.xticks(x, mods)
    xticks(rotation=90)
    ax[:set_xlim]([0, length(mods)+1])

    tight_layout()
    subplots_adjust(hspace=0.15) 
end

## CODES

0) obs.jl
    observables as a function of euclidean time x0 and model average. Hardcoded to explore two fit functions in model av: constant (plateau) and exponential + constant. You need to give two set of lower bounds and upper bounds fit interval (one for exponential and other for plateau), e.g. tm=[[1,2,3], [1,2,3]], tM=[[100,101,102], [100,101,102]]. First comes the exponential + constant fit, then the plateau fit

1) wil_sym.jl 
    computes (t0, mpi, mk=mpi, m12, m13=m12, fpi, fk=fpi) in wilson regularization for symmetric point ensembles
    improves and renormalizes (ZA) fpi & fk
    does not improve nor renormalize m12 & m13
    stores in bdio (t0, mpi, mk=mpi, m12, m13=m12, fpi, fk=fpi)

2) wil_constTr.jl 
    computes (t0, mpi, mk, m12, m13, fpi, fk) in wilson regularization for non symmetric point ensembles
    improves and renormalizes (ZA) fpi & fk
    does not improve nor renormalize m12 & m13
    stores in bdio (t0, mpi, mk, m12, m13, fpi, fk)

3) to mass shift, read der_1q.bdio & der_sea_1q -> derivatives of sqrt(t0)fpi & sqrt(t0)fk are computed with improvement and renormalization
                                                   (ZA wilson, 1 wtm)
                                                   derivatives of sqrt(t0)m12 & sqrt(t0)m13 are computed without improvement but with renormalization (ZA/ZP wilson, 1/ZP wtm)

4) matching_sym.jl & matching_constTr.jl
    read unshifted wilson & wtm data
        when reading, fpi & fk are renormalized and improved, m12 & m13 are unrenormalized and unimproved
    mass shift wilson & wtm
    matching & full twist
    interpolate decay constants to matching & full twist point
    stores in bdio wilson (t0_sh, phi2_sh, fpik_sh), phi4_sh can be known from in.jl, fpik_sh is renormalized and improved and in units of sqrt(t0)
                   wtm matched (kappa*, mul*, muls*, fpik*), fpik* in units of sqrt(t0)
    ## check improvement and renormalization of pcac masses when storing _sh in bdio

## BDIOs

1) ens_obs_wil_un.bdio
    saved (t0,mpi,mk,m12,m13,fpi,fk) wilson, fpi & fk improved and renormalized, m12 & m13 unimproved and unrenormalized

2) ens_"obs"_tm_un.bdio
    saved the grid values of "obs"=(t0,mpi,mk,m12,m13,fpi,fk) wtm
        wtm -> no improvement in fpi, fk, m12, m13
        wtm -> no renormalization in fpi & fk
        saved m12 & m13 are bare pcac quark masses, unrenormalized

3) ens_obs_wil_sh_phi4=...bdio
    saved (t0, phi2, sqrt(t0)fpik) wilson shifted to phi4 in name of the file

4) ens_obs_tm_sh_phi4=...bdio
    saved (kappa*, mul*, muls*, sqrt(t0)fpik) wtm matched (*) shifted to phi4 in name of the file

5) der_1q.bdio
    saved parameters needed for mass shifting different observables
    parameters
        1:3   -> sqrt(t0)fpik wilson (improved and renormalized)
        4:6   -> phi2 wilson
        7:9   -> t0
        10:12 -> sqrt(t0)fpik wtm (no need for improvement or renormalization since wtm)
        13:15 -> phi2 wtm (grid)
        16:18 -> phi4 wtm (grid)
        19:23 -> sqrt(t0)m12pcac/ZP wtm (grid)
    ## update w/ m12 wilson, m13 wilson, sqrt(t0)fpi wilson & wtm, sqrt(t0)fk wilson & wtm

6) der_ens_1q.bdio
    saved the mass derivatives of wilson observables: (valence + sea) derivative of observable wrt phi4_wilson
    parameters
        1 -> derivative of sqrt(t0)fpik wilson (improved and renormalized) wrt phi4_wilson
        2 -> ... of phi2 wilson ...
        3 -> t0
        4 -> sqrt(t0)fpi wilson (impr and ren)
        5 -> sqrt(t0)fk wilson (impr and ren)
        6 -> sqrt(t0)m12pcac wilson -> unimproved and unrenormalized -> before extrapolating in beta and phi2, need to renormalize (ZA/ZP) to be phys.
        7 -> sqrt(t0)m13pcac wilson ...

7) der_sea_ens_1q.bdio
    mass derivative of wtm observables: (sea) derivative wrt phi4_wilson -> I take one point in the grid for each ensemble, since phi2~phi2_wilson
    parameters
        1 -> sqrt(t0)fpik wtm (wtm thus unimpr and unren, but no need, already physical)
        2 -> phi2 wtm
        3 -> phi4 wtm
        4 -> sqrt(t0)m12pcac wtm -> before extrapolating in beta and phi2, need to renormalize with 1/ZP to be physical
        5 -> sqrt(t0)fpi wtm
        6 -> sqrt(t0)fk wtm

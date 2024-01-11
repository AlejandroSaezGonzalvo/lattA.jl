phi2 = [8 * t0 * mpi[i] ^ 2 for i in 1:length(mpi)]
phi4 = Array{uwreal,1}()
fpik = Array{uwreal,1}()
c=0
for i in 1:length(mpi)
    for j in 2*i-1+c:2*i+1+c
        push!(phi4, 8 * t0 * (mk[j] ^ 2 + 0.5 * mpi[i] ^ 2))
        push!(fpik, sqrt(t0) * 2/3 * (fk[j] + 0.5 * fpi[i]))
    end
    c+=1
end
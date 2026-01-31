module ramp
push!(LOAD_PATH, pwd())
using LinearAlgebra
using QuantumOptics
using DifferentialEquations
include("diagonalization.jl")
using .diagonalization
export adiabaticramp


function adiabaticramp(psi0,psif,HH,K,JJ,epsilon,tau,Ntau,acc)
    tmax = Ntau*tau
    psi0t = transpose(conj(psi0))
    psift = transpose(conj(psif))
    times=(0.0,tmax)
    ev = eigvals(HH)
    tint=tau/acc
    psit    = psi0
    psitold = psi0
    open("output/fidelities_ramp.dat","w") do io
    for k in 0:Ntau
        for ki in 1:100
        psitold = psit    
        psit = evolt(psitold,HH,tint,K,epsilon,tau,tmax,k*tau + (ki-1)*tint)
        psit =  ComplexF64.(psit)
        end
        ovl1 = psi0t*psit
        ovl2 = psift*psit
        ovl3 = psi0t*psif
        println(io,k*tau," ",abs2(ovl1)," ",abs2(ovl2)," ",abs2(ovl3))
    end
    end
  return psit
end

function evolt(psi0,HH,dt,K,epsilon,tau,tramp,t)
    U = exp(-im*dt*(HH + (epsilon/2)*(1-cos(pi*(t+dt/2)/tramp))*K*cos(2*pi*(t+dt/2)/tau)))
    psit = U*psi0
    return psit
end


end

module qfi
push!(LOAD_PATH, pwd())
using LinearAlgebra
using QuantumOptics
include("diagonalization.jl")
include("wigner.jl")
using .diagonalization
using .wigner

#println("hereee")

function qfirotationJz(psi,J)
    Jz = diagonalization.matrixJz(J) 
    psit = psi'
    jz  = real(psit*Jz*psi)
    jz2 = real(psit*(Jz^2)*psi)
    varjz = jz2 - jz^2
    qfi = 4*varjz
    nqfi = qfi/(2*J-2*jz^2/J)
    return nqfi
end

function qfirotationJy(psi,J)
    Jp = diagonalization.matrixJp(J)
    Jm = Jp'
    Jy = (1/(2*im))*(Jp-Jm)
    psit = psi'
    jy  = real(psit*Jy*psi)
    jy2 = real(psit*(Jy^2)*psi)
    varjy = jy2 - jy^2
    qfi = 4*varjy
    nqfi = qfi/(2*J-2*jy^2/J)
    return nqfi
end


function qfi_mixed(rho,A,J)
  rhodiag = eigen(rho)
  lam = rhodiag.values
  rhovectors = rhodiag.vectors
  qfi2t = 0 + 0*im
  qfi1t = 0 + 0*im
  epsilon=0.0000001    
  for k in 1:length(lam)
      for l in 1:k-1
       kstate = [rhovectors[i,k] for i in 1:length(lam)]
       lstate = [rhovectors[i,l] for i in 1:length(lam)]
       lstatead = transpose(conj(lstate))
       sd = lstatead*A*kstate
       if (lam[k]+lam[l])>epsilon
         qfiinst2t = 8*((lam[l]*lam[k])/(lam[k]+lam[l]))*abs2(sd) 
         qfi2t = qfi2t + qfiinst2t
       end
      end
      for l in k+1:length(lam)
       kstate = [rhovectors[i,k] for i in 1:length(lam)]
       lstate = [rhovectors[i,l] for i in 1:length(lam)]
       lstatead = transpose(conj(lstate))
       sd = lstatead*A*kstate
       if (lam[k]+lam[l])>epsilon
         qfiinst2t = 8*((lam[l]*lam[k])/(lam[k]+lam[l]))*abs2(sd) 
         qfi2t = qfi2t + qfiinst2t
       end
     end
     kstate = [rhovectors[i,k] for i in 1:length(lam)]
     kstatead = transpose(conj(kstate))
     sq2ev = kstatead*(A^2)*kstate
     sqev2 = abs2(kstatead*(A)*kstate)
     qfiinst1t = 4*lam[k]*(sq2ev -sqev2)
     qfi1t = qfi1t + qfiinst1t
  end
  jexp = tr(rho*A)
  return real(qfi1t - 1.0*qfi2t)/(2*J-2*jexp^2/J)
end



# ----- Imput parameters  ----
#J=20                       # System size
#NN=100                     # Size of the Grid for phase space
#name1="../output/wignertest_kicked.dat"     # Wigner output file
#alpha = 0.0 + 0.5*im
#theta = 1.0
# ------------------------- ----

#N=2*J+1
#psi0 = zeros(ComplexF64, N)
#psi0[1]= 1.0 + 0.0*im

#Jp = diagonalization.matrixJp(J)
#Jz = diagonalization.matrixJz(J)
#Jx = diagonalization.matrixJx(J)
#Jm = Jp'
#Jy = (1/(2*im))*(Jp-Jm)
#psic = exp(alpha*Jp - conj(alpha)*Jm)*psi0
#psicat = (1/2^(1/2))*(exp(alpha*Jp - conj(alpha)*Jm)*psi0 + exp((-alpha)*Jp - conj(-(alpha))*Jm)*psi0)
#psic1 = exp(alpha*Jp - conj(alpha)*Jm)*psi0
#psic2 = exp(-alpha*Jp - conj(-alpha)*Jm)*psi0
#psi = exp(-im*Jy*theta)*psi
#psiqo = wigner.buildingstate(psi,J)
#qfiy = qfirotationJy(psic,J)
#println("qfi: ", qfiy)
#psicatt = psicat'
#rhocat = psicat*psicatt
#qfiy_cat = qfi_mixed(rhocat,Jy,J)
#println("qfi_cat : ",qfiy_cat)
#psic1t = psic1'
#psic2t = psic2'
#rhoc1 = psic1*psic1t
#rhoc2 = psic2*psic2t
#rho_mix = (1/2)*rhoc1 + (1/2)*rhoc2
#qfiy_mix = qfi_mixed(rho_mix,Jy,J)
#println("qfi_cat : ",qfiy_mix)
#rhoqo_mix = wigner.buildingstaterho(rho_mix,J)
#psiqo_cat = wigner.buildingstate(psicat,J)
#println("Normalized QFI(Jz): ",qfir)
#wtest = wigner.wignerf(rhoqo_mix,NN,name1)
#wtest = wigner.wignerf(psiqo_cat,NN,name1)

#println("done")



end

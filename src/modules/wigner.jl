module wigner
push!(LOAD_PATH, pwd())
using LinearAlgebra
using QuantumOptics
include("diagonalization.jl")
include("troterization.jl")
using .diagonalization
using .troterization
export husimif
export wignerf
export buildingstate


function husimif(psi,JJ,NN,name)
  thetalist = [i*(pi/NN) for i in 1:NN]
  philist = [i*(2*pi/NN) for i in 1:NN]
  psiad = conj(transpose(psi))
  rho = psi*psiad
  Jpop = diagonalization.matrixJp(JJ)
  mJstate =[0.0 for i in -JJ:JJ]
  mJstate[1] = 1.0
  open(name,"w") do io
  for theta in thetalist
    for phi in philist
      al = tan(theta/2)*exp(-im*phi)
      alket = (1/(1+abs2(al))^(JJ))*exp(al*Jpop)*mJstate
      albra= conj(transpose(alket))
      hus= (1/pi)*albra*rho*alket
      println(io,theta," ",phi," ",real(hus[1]))
    end
  end
  end
  return "done"
end

function wignerf(qostate,NN,name)
  thetalist = [i*(pi/NN) for i in 1:NN]
  philist = [i*(2*pi/NN) for i in 1:NN]
  open(name,"w") do io
  for i in thetalist
    for j in philist
      wig = wignersu2(qostate,i,j)
      println(io,i," ",j," ",real(wig))
    end
  end
  end
  return "done"
end


function buildingstate(lista,jj)
  ba=SpinBasis(jj)
  psiinst = spindown(ba)
  psif = lista[1]*psiinst
  for i in 2:length(lista)
    psiinst=sigmap(ba)*psiinst/(jj*(jj+1)-((-jj+i-2)*((-jj+i-2)+1)))^(1/2)
    psif = psif + lista[i]*psiinst
  end
  return psif
end

function HLMG_qop(J,ep,gx,gy)
  ba=SpinBasis(J)
  v=ep*(gx-gy)/(2.0*(2.0*J-1.0))
  w=ep*(gx+gy)/(2.0*(2.0*J-1.0))
  H = (ep/2)*sigmaz(ba) + 1.0*(v/2)*((sigmap(ba))^2+(sigmam(ba))^2) + 1.0*(w/2)*(sigmap(ba)*sigmam(ba) + sigmam(ba)*sigmap(ba))
  return H
end


#J=50
#ep=1.0
#gx=-3.0
#gy=-9.0
#NN=100
#name1="wignertest.dat"
#name2="husimitest.dat"
#k = 3 # state of interest

#HH0 = diagonalization.matrixH00(J,ep,gx,gy)
#ev0 = eigvals(HH0)
#println("flag 1: ",ev0[k])
# building Hamiltonian Quantum optics
#Hqo = HLMG_qop(J,ep,gx,gy)
#ev2 = eigenstates(dense(Hqo))
#println("flag 2: ",ev2[1][k])
#
#evecs = eigvecs(HH0)
#listags = [evecs[i,k] for i in 1:(2*J+1)]
#psi0 = buildingstate(listags,J)
#testev = expect(Hqo,psi0)
#println("flag 3: ",real(testev))
#htest = husimif(psi0,NN,name2)
#println("Husimi done")
#wtest = wignerf(psi0,NN,name1)



end

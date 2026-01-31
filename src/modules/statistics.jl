module statistics
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("diagonalization.jl")
include("troterization.jl")
using .diagonalization
using .troterization



function analysisH(N,om,r,lambda,delta,nn,nu,chi)
   eigvs=diagonalization.diagonalize(N,om,r,lambda,delta)
   open("DensityOfStates_output.dat","w") do ioa
   for i in 1:trunc(Int64,length(eigvs[1])/2-1)
      println(ioa,(eigvs[1][2*i+1]+eigvs[1][2*i])/2," ",1.0/(eigvs[1][2*i+1]-eigvs[1][2*i]))
   end
   end
   println("See file DensityOfStates_output.dat for Density of states (DOS) ")
   floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu)
   eigvecsf=eigvecs(floquet)
   ham=diagonalization.hamiltonian(N,om,r,lambda,delta)
   evfvec = zeros(0)
   open("levels_output.dat","w") do io
   for i in 1:2*(N+1)
     eigvecf=[eigvecsf[j,i] for j in 1:2*(N+1)]
     eigvecf_t=Array{Complex{Float64}}(undef,1,length(eigvecf)) 
     for k in 1:length(eigvecf)
       eigvecf_t[1,k]=conj(eigvecf[k])
     end
     evf=eigvecf_t*ham*eigvecf
     append!(evfvec, real(evf) )
   end
   evfvecs=sort(evfvec)
   for i in 1:length(evfvec)
     println(io,i," ",eigvs[1][i]," ",evfvecs[i])
   end
   end
   return "Done"
end


function parameter_r(J,ep,gx,gy,nn,nu,chi)
   T=2*pi/nu
   #floquet=troterization.troter(N,nn,r,om,lambda,delta,chi,nu)
   H0=diagonalization.matrixH0(J,ep,gx,gy)
   Jzdiag=[(-J+2*(jj)) for jj in 0:J]
   Jz=Array(Diagonal(Jzdiag))
   floquet=exp(-im*chi*Jz)*exp(-im*H0*T)
   floquet2=troterization.troter(H0,Jz,J,nn,chi,nu)
   eigvalsf=eigvals(floquet)
   eigvalsf2=eigvals(floquet2)
   qs=zeros(0)
   qs2=zeros(0)
   for i in 1:length(eigvalsf)
      f1p=mod(real(log(eigvalsf[i])/(-im)),2*pi)
      append!(qs,f1p)
      #println(f1p)
      f2p=mod(real(log(eigvalsf2[i])/(-im)),2*pi)
      append!(qs2, f2p )
   end
   qs=sort(qs)
   qs2=sort(qs2)
   #println(qs)
   sn=[qs[n+1]-qs[n] for n in 1:(length(qs)-1)]
   sn2=[qs2[n+1]-qs2[n] for n in 1:(length(qs2)-1)]
   #sn=[qs[2*n+1]-qs[2*n-1] for n in 1:trunc(Int64,length(qs)/2-1)]
   #
   fac=sum(sn)/length(sn)
   fac2=sum(sn2)/length(sn2)
   sn=(1/fac)*(sn)
   sn2=(1/fac2)*(sn2)
   #println(sn)
   rn=zeros(0)
   rn2=zeros(0)
   #open("quasienergies_spacing.dat","w") do io
   for n in 1:length(sn)-1
       if sn[n]>sn[n+1]
         append!(rn, sn[n+1]/sn[n])
       end
       if sn[n]<sn[n+1]
         append!(rn, sn[n]/sn[n+1])
       end
   #  println(io," ",sn[n]) 
   end
   for n in 1:length(sn2)-1
       if sn2[n]>sn2[n+1]
         append!(rn2, sn2[n+1]/sn2[n])
       end
       if sn2[n]<sn2[n+1]
         append!(rn2, sn2[n]/sn2[n+1])
       end
#     println(io," ",sn[n]) 
   end
#   end
   rpar=sum(rn)/length(rn)
   rpar2=sum(rn2)/length(rn2)
   return [rpar,rpar2]
   end

function sortingvecf(floquet,H0,JJ)
   eigvecsf=eigvecs(floquet)
   vals = eigvals(H0)
   evfvec = zeros(0)
   for i in 1:length(vals)
     eigvecf=[eigvecsf[j,i] for j in 1:length(vals)]
     eigvecf_t=Array{Complex{Float64}}(undef,1,length(eigvecf)) 
     for k in 1:length(eigvecf)
       eigvecf_t[1,k]=conj(eigvecf[k])
     end
     evf=eigvecf_t*H0*eigvecf
     append!(evfvec, real(evf) )
   end
   evford=sortperm(evfvec)
   return evford
end

function expectation(floquet,H0,J,ep,gx,gy,chi,T,k)
   #H0=diagonalization.matrixH0(J,ep,gx,gy)
   ev0=eigvals(H0)
   #Jzdiag=[(-J+2*(jj)) for jj in 0:J]
   Jxop=diagonalization.matrixJx(J)
   #floquet=exp(-im*chi*Jz)*exp(-im*H0*T)
   eigvecsfloquet=eigvecs(floquet)
   so=sortingvecf(floquet,H0,J)
   listvec=[eigvecsfloquet[i,so[k]] for i in 1:length(ev0)]
   listvect=conj(transpose(listvec))
   evf = listvect*H0*listvec
   qfiJx = 4*(listvect*(Jxop^2)*listvec - (listvect*(Jxop)*listvec)^2) 
   evalfloquet = real(evf)
   nqfiJx = qfiJx /(2*J+1) 
   return [ev0[k],evalfloquet,nqfiJx]
   end

#intchi=0.02
#chilist=[i*intchi for i in 1:500]
#for cc in chilist
#   rtest=parameter_r(J,ep,gx,gy,N,nu,cc)
#   println(cc," ",rtest[1]," ", rtest[2])
#end


#J=50
#ep=1.0
#gx=-0.95
#gy=3*gx
#nu=1
#N=20
#chi=0.001

#open("resonances.dat","w") do io

#T=0.1
#tint = 0.01
#Nmax=400
#for ik in 1:Nmax
#   test1 = expectation(J,ep,gx,gy,chi,T,1)
#   test2 = expectation(J,ep,gx,gy,chi,T,2)
#   test3 = expectation(J,ep,gx,gy,chi,T,3)
#   test4 = expectation(J,ep,gx,gy,chi,T,4)
#   test5 = expectation(J,ep,gx,gy,chi,T,5)
#   println(io,T," ",test1[2]," ",test2[2]," ",test3[2]," ",test4[2]," ",test5[2])
#   println(io,T," ",test1[2])
#   println(ik/Nmax)
#   T=T+tint
#end

#end

end

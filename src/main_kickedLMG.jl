module main_kickedLMG
push!(LOAD_PATH, pwd())
using LinearAlgebra
using QuantumOptics
include("modules/diagonalization.jl")
include("modules/troterization.jl")
include("modules/wigner.jl")
include("modules/statistics.jl")
using .diagonalization
using .troterization
using .wigner
using .statistics




# ----- Imput parameters  ----
k = 1                      # state of interest
J=50                       # System size
ep=1.0                     # LMG parameter
gx=-0.95                   # LMG parameter
#gx=-0.25                   # LMG parameter
gy= 3*gx                   # LMG parameter
epsilon = 0.01             # strength of the Kick
#tau = 2.89134             # period of the kick
tau = 2.89134/2            # period of the kick
NN=100                     # Size of the Grid
name1="output/wignertest_kicked.dat"     # Wigner output file
name2="output/husimitest_kicked.dat"     # Husimi output file
# ------------------------- ----


#---------------- Main body -----------------------
#--------------------------------------------------
#---------  Building Floquet operator -------------#
HH0 = diagonalization.matrixH00(J,ep,gx,gy)
ev0 = eigvals(HH0)
states0 = eigvecs(HH0)
vec0=[states0[i,k] for i in 1:length(ev0)]
vec0t=conj(transpose(vec0))
Jz =  diagonalization.matrixJz(J)
Jx =  diagonalization.matrixJx(J)
Kop = Jz+Jx
Floquet = exp(-im*tau*HH0)*exp(-im*epsilon*Kop)
fstates = eigvecs(Floquet)
fev     = eigvals(Floquet)
println("-> Floquet operator obtained")
#---   sorting Floquet states   -------------#
sortf   = statistics.sortingvecf(Floquet,HH0,J)
listvec=[fstates[i,sortf[k]] for i in 1:length(fev)]
println("-> Stationary state ",k,"-th obtained")
# building floquet state for QuantumOptics library
psi = wigner.buildingstate(listvec,J)
# calculating Husimi and Wigner using Quantum optics library
htest = wigner.husimif(listvec,J,NN,name2)
println("-> Husimi function obtained")
wtest = wigner.wignerf(psi,NN,name1)
println("-> Wigner function obtained")
# calculating <F_k|H0|F_k>
psif = listvec
psift = conj(transpose(psif))
fexpval = real(psift*HH0*psif)
qfi_F = 4*real(psift*(Jz^2)*psif-(psift*(Jz)*psif)^2)/J
qfi_E = 4*real(vec0t*(Jz^2)*vec0-(vec0t*(Jz)*vec0)^2)/J
#------------------------------------------------------


#
open("output/resonances_kicked.dat","w") do io
tint=0.01
T = [i*tint for i in 20:500]
for tinst in T
   Floquet = exp(-im*tinst*HH0)*exp(-im*epsilon*Kop)
   test1 = statistics.expectation(Floquet,HH0,J,ep,gx,gy,epsilon,tinst,k)
   println(io,tinst," ",real(test1[2])/J," ",real(test1[3])/J,)
   println(tinst/T[length(T)])
end
end


println("See file output/resonances_kicked.dat and see: - First column:     Tau")
println("                                        - Second column: <F_k|H0|F_k>/J")
println("                                        - Third column:  QFI[rho_k,Jz]/J")

# ----  Printing results ------------
println("----- Results ----------")
println("Floquet stationary state k: ",k)
println("<E_k|H0|E_k>/J : ",ev0[k]/J)
println("QFI[|E_k><E_k|,Jz]/J: ",qfi_E/J)
println("<F_k|H0|F_k>/J : ",fexpval/J)
println("QFI[|F_k><F_k|,Jz]/J: ",qfi_F/J)
println("Go to file ",name2," for Husimi function")
println("Go to file ",name1," for Wigner function")
println("-----------------------" )
#-------------------------------------

end

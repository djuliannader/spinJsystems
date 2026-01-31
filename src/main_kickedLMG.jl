module main_kickedLMG
push!(LOAD_PATH, pwd())
using LinearAlgebra
using QuantumOptics
include("modules/diagonalization.jl")
include("modules/troterization.jl")
include("modules/wigner.jl")
include("modules/statistics.jl")
include("modules/qfi.jl")
using .diagonalization
using .troterization
using .wigner
using .statistics
using .qfi




# ----- Imput parameters  ----
k = 1                      # state of interest
J=50                       # System size
ep=1.0                     # LMG parameter
gx= -0.6                   # LMG parameter
gy= 0                      # LMG parameter
epsilon = 0.01             # strength of the Kick
#tau = 9.884740906201314  # period of the kick
tau = 9.88
NN=100                     # Size of the Grid for phase space
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
println("-> Floquet state ",k,"-th obtained")
# building Floquet state for QuantumOptics library
psi = wigner.buildingstate(listvec,J)
# calculating Husimi and Wigner using Quantum optics library
psiH0 = listvec
psiH0t = psiH0'
rho_H0 = psiH0*psiH0t
htest = wigner.husimif(rho_H0,J,NN,name2)
println("-> Husimi function obtained")
wtest = wigner.wignerf(psi,NN,name1)
println("-> Wigner function obtained")
# calculating <F_k|H0|F_k>
psif = listvec
psift = conj(transpose(psif))
fexpval = real(psift*HH0*psif)
#qfi_F = 4*real(psift*(Jz^2)*psif-(psift*(Jz)*psif)^2)/J
#qfi_E = 4*real(vec0t*(Jz^2)*vec0-(vec0t*(Jz)*vec0)^2)/J
qfi_F = qfi.qfirotationJz(psif,J)
qfi_E = qfi.qfirotationJz(vec0,J)
qfi_F2 = qfi.qfirotationJy(psif,J)
qfi_E2 = qfi.qfirotationJy(vec0,J)
#------------------------------------------------------




#
open("output/resonances_kicked.dat","w") do io
tint=0.01
T = [i*tint for i in 20:1500]
for tinst in T
   Floquet = exp(-im*tinst*HH0)*exp(-im*epsilon*Kop)
   test1 = statistics.expectation(Floquet,HH0,J,ep,gx,gy,epsilon,tinst,k)
   sortf   = statistics.sortingvecf(Floquet,HH0,J)
   fstates = eigvecs(Floquet) 
   listvectau=[fstates[i,sortf[k]] for i in 1:length(fev)]
   qfi_tauz = qfi.qfirotationJz(listvectau,J)
   qfi_tauy = qfi.qfirotationJy(listvectau,J) 
   println(io,tinst," ",real(test1[2])/J," ",qfi_tauz," ",qfi_tauy)
   println(tinst/T[length(T)])
end
end

# ----  Printing results ------------
println("----- Results ----------")
println("Floquet stationary state k: ",k)
println("<E_k|H0|E_k>/J : ",ev0[k]/J)
println("Normalized QFI[|E_k><E_k|,Jz]/(2J-2<Jz>/J): ",qfi_E)
println("Normalized QFI[|E_k><E_k|,Jy]/(2J-2<Jz>/J): ",qfi_E2)
println("<F_k|H0|F_k>/J : ",fexpval/J)
println("Normalized QFI[|F_k><F_k|,Jz]/(2J-2<Jz>/J): ",qfi_F)
println("Normalized QFI[|F_k><F_k|,Jy]/(2J-2<Jz>/J): ",qfi_F2)
println("Go to file ",name2," for Husimi function")
println("Go to file ",name1," for Wigner function")
println("-----------------------" )
#-------------------------------------


println("See file output/resonances_kicked.dat and see: - First column:     Tau")
println("                                        - Second column: <F_k|H0|F_k>/J")
println("                                        - Third column:  QFI[rho_k,Jx]/(2J+1)")



end

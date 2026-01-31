module main_adiabaticramp
push!(LOAD_PATH, pwd())
using LinearAlgebra
using QuantumOptics
include("modules/diagonalization.jl")
include("modules/troterization.jl")
include("modules/wigner.jl")
include("modules/statistics.jl")
include("modules/ramp.jl")
using .diagonalization
using .troterization
using .wigner
using .statistics
using .ramp




# ----- Imput parameters  ----
k = 1                      # state of interest
J=50                       # System size
ep=1.0                     # LMG parameter
gx=-0.95                   # LMG parameter
gy= 3*gx                   # LMG parameter
epsilon = 0.01             # amplitude of the drive
#tau = 2.89134              # period of the drive
tau = 2.5                 # period of the drive
Nramp = 30                 # periods of the ramp
NN=100                     # Size of the Grid for phase space
nfloquet=200               # Number of time intervals for the Trotter-Suzuki decomposition
Ndt = 100                  # Number of subintervals of the modulation
name1="output/wignertest_ramp.dat"     # Wigner output file
name2="output/husimitest_ramp.dat"     # Wigner output file
# ------------------------- ----


#---------------- Main body -----------------------
#--------------------------------------------------
#---------  Building initial state -------------#
HH0 = diagonalization.matrixH00(J,ep,gx,gy)
ev0 = eigvals(HH0)
states0 = eigvecs(HH0)
psi0=[states0[i,k]+0.0*im for i in 1:length(ev0)]
# - Building Floquet state
Jz =  diagonalization.matrixJz(J)
Jx =  diagonalization.matrixJx(J)
Kop = Jz+Jx
Floquet = troterization.troter(HH0,Kop,J,nfloquet,epsilon,tau)
fstates = eigvecs(Floquet)
fev     = eigvals(Floquet)
println("-> Floquet operator obtained")
# Building target (Floquet) state
sortf   = statistics.sortingvecf(Floquet,HH0,J)
psif=[fstates[i,sortf[k]] for i in 1:length(fev)]
println("-> Target Floquet state ",k,"-th obtained")
# Adiabatic ramp
psit = ramp.adiabaticramp(psi0,psif,HH0,Kop,J,epsilon,tau,Nramp,Ndt)
# building final state for QuantumOptics library
psit_qo = wigner.buildingstate(psit,J)
# calculating Husimi and Wigner using Quantum optics library
htest = wigner.husimif(psit,J,NN,name2)
println("-> Husimi function obtained")
wtest = wigner.wignerf(psit_qo,NN,name1)
println("-> Wigner function obtained")
# calculating QFI
#psitt = conj(transpose(psit))
#fexpval = real(psitt*HH0*psit)
#qfi_F = 4*real(psitt*(Jz^2)*psit-(psitt*(Jz)*psit)^2)/J
#------------------------------------------------------



# ----  Printing results ------------
println("----- Results ----------")
println("Go to file ",name2," for Husimi function")
println("Go to file ",name1," for Wigner function")
println("-----------------------" )
#-------------------------------------


println("See file output/fidelities_ramp.dat and see: - First column:     time")
println("                                        - Second column: with respect to the initial state <E_k|psi(t)> ")
println("                                        - Third column:  with respect to the target state  <F_k|psi(f)> ")



end

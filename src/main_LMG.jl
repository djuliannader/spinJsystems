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
gx=0.75                    # LMG parameter
gy=0.0                    # LMG parameter
NN=100                     # Size of the Grid
name1="output/wignertest.dat"     # Wigner output file
name2="output/husimitest.dat"     # Husimi output file
# ------------------------- ----


#---------------- Main body -----------------------
#--------------------------------------------------
#---------  Building Floquet operator -------------#
HH0 = diagonalization.matrixH00(J,ep,gx,gy)
ev0 = eigvals(HH0)
states = eigvecs(HH0)
println("-> Eigenstates obtained")
# building state for QuantumOptics library
liststate = [states[i,k] for i in 1:(2*J+1)]
psi = wigner.buildingstate(liststate,J)
# calculating Husimi and Wigner using Quantum optics library
htest = wigner.husimif(psi,NN,name2)
println("-> Husimi function obtained")
wtest = wigner.wignerf(psi,NN,name1)
println("-> Wigner function obtained")
#------------------------------------------------------



# ----  Printing results ------------
println("----- Results ----------")
println("LMG stationary state: ",k)
println("<E_k|H0|E_k>/J: ",ev0[k]/J)
println("Go to file ",name2," for Husimi function")
println("Go to file ",name1," for Wigner function")
println("-----------------------" )
#-------------------------------------




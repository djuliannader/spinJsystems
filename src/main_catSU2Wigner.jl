module qfi
push!(LOAD_PATH, pwd())
using LinearAlgebra
using QuantumOptics
include("modules/diagonalization.jl")
include("modules/wigner.jl")
include("modules/qfi.jl")
using .diagonalization
using .wigner
using .qfi




# ----- Imput parameters  ----
J=20                                         # System size
NN=100                                       # Size of the Grid for phase space
alpha = 0.3 + 0.0*im                         # displacement parameter
flag = 0                                     # (1) Wigner function (0) not Wigner function
name1="output/wignertest_cats.dat"         # Wigner output file 1
name2="output/wignertest_cata.dat"         # Wigner output file 2
name3="output/wignertest_mixs.dat"         # Wigner output file 3
name4="output/wignertest_mixa.dat"         # Wigner output file 4
name5="output/husimitest_cats.dat"         # Wigner output file 5
name6="output/husimitest_cata.dat"         # Wigner output file 6
name7="output/husimitest_mixs.dat"         # Wigner output file 7
name8="output/husimitest_mixa.dat"         # Wigner output file 8
# ------------------------- ----

# Body
# ---------------------------------------------
# generating state of minimal projection |J,-J>
N=2*J+1
psi0 = zeros(ComplexF64, N)
psi0[1]= 1.0 + 0.0*im
#
# generating necessary operators
Jp = diagonalization.matrixJp(J)
Jz = diagonalization.matrixJz(J)
Jm = Jp'
Jx = diagonalization.matrixJx(J)
Jy = (1/(2*im))*(Jp-Jm)
#
# building cat and mixed states
psicats = (1/2^(1/2))*(exp(alpha*Jp - conj(alpha)*Jm)*psi0 + exp((-alpha)*Jp - conj(-(alpha))*Jm)*psi0)
psicatst = psicats'
rho_cats = psicats*psicatst
rhoqo_cats = wigner.buildingstaterho(rho_cats,J)
psicata = (1/2^(1/2))*(psi0 + exp((2^(1/2)*alpha)*Jp - conj((2^(1/2)*alpha))*Jm)*psi0)
psicatat = psicata'
rho_cata = psicata*psicatat
rhoqo_cata = wigner.buildingstaterho(rho_cata,J)
psic1 = exp(alpha*Jp - conj(alpha)*Jm)*psi0
psic2 = exp(-alpha*Jp - conj(-alpha)*Jm)*psi0
psic3 = exp((2^(1/2)*alpha)*Jp - conj(2^(1/2)*alpha)*Jm)*psi0
psic1t = psic1'
psic2t = psic2'
psic3t = psic3'
psi0t  = psi0'
rho_c1 = psic1*psic1t
rho_c2 = psic2*psic2t
rho_c3 = psic3*psic3t
rho_0  = psi0*psi0t
rho_mixs = (1/2)*rho_c1 + (1/2)*rho_c2
rho_mixa = (1/2)*rho_0 + (1/2)*rho_c3
rhoqo_mixs = wigner.buildingstaterho(rho_mixs,J)
rhoqo_mixa = wigner.buildingstaterho(rho_mixa,J)
# calculatin Expectation jz
jz_cats = tr(rho_cats*Jz)
jz_cata = tr(rho_cata*Jz)
jz_mixs = tr(rho_mixs*Jz)
jz_mixa = tr(rho_mixa*Jz)
# calculating QFI
qfi_cats_Jz = qfi.qfi_mixed(rho_cats,Jx,J)
qfi_cata_Jz = qfi.qfi_mixed(rho_cata,Jx,J)
qfi_mixs_Jz = qfi.qfi_mixed(rho_mixs,Jx,J)
qfi_mixa_Jz = qfi.qfi_mixed(rho_mixa,Jx,J)



# printing results

println("---- Expectation values <J_z> ---")
println("  Cat state     (s) : ",jz_cats)
println("  Cat state     (a) : ",jz_cata)
println("  Mixture state (s) : ",jz_mixs)
println("  Mixture state (a) : ",jz_mixa)
println("----    QFI(J_z)      ---")
println("  Cat state     (s) : ",real(qfi_cats_Jz))
println("  Cat state     (a) : ",real(qfi_cata_Jz))
println("  Mixture state (s) : ",real(qfi_mixs_Jz))
println("  Mixture state (a) : ",real(qfi_mixa_Jz))




if flag == 1
    wtest1 = wigner.wignerf(rhoqo_cats,NN,name1)
    wtest2 = wigner.wignerf(rhoqo_cata,NN,name2)
    wtest3 = wigner.wignerf(rhoqo_mixs,NN,name3)
    wtest4 = wigner.wignerf(rhoqo_mixa,NN,name4)
end

    htest1 = wigner.husimif(rho_cats,J,NN,name5)
    htest2 = wigner.husimif(rho_cata,J,NN,name6)
    htest3 = wigner.husimif(rho_mixs,J,NN,name7)
    htest4 = wigner.husimif(rho_mixa,J,NN,name8)



end

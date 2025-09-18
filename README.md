# driven LMG
Numerical study of the Floquet states of the Lipkin Meshkov Glick Model with a periodic drive

clone as:

```bash
git clone https://github.com/djuliannader/drivenLMG.git
```

## Data 

This repository provides data and codes for the papers for the Floquet states of the LMG by J. Segura, et al.
The paper are available on PRA.
The repository is currently being developed; please await its completion1.  The files contain a Mathematica notebook, necessary for the analysis of the results.


## Environment Requirements  

To run this repository, please make sure the following environment is available:

- LinearAlgebra (stdlib)  
- DifferentialEquations
- HCubature
- QuantumOptics

## Usage

- For the stationary states of the time-independent LMG manipulate parameters, move to the file /src/main_LMG.jl manipulate parameters and run as:

```bash
julia main_LMG.jl
```

Cite as: D. Nader, C. Gonzalez-Rodrıguez, and S. Lerma-Hernandez, Avoided crossings and dynamical tunneling
close to excited-state quantum phase transitions, Physical Review E. (2021).

- For the Floquet states of the kicked LMG manipulate, move to the file /src/main_kickedLMG.jl manipulate parameters and run as:

```bash
julia main_kickedLMG.jl
```

Cite as: J. Segura-Landa et al, Quantum precursors to Kolmogorov-Arnold-Moser theorem in Floquet spin-J systems, arXiv:2504.13257 (2025)


- For the Floquet states of the driven LMG manipulate, move to the directory /src/main_drivenLMG.jl manipulate parameters and run as:

```bash
julia main_drivenLMG.jl
```

Cite as: (In progress)



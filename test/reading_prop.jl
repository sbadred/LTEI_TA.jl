using LTEI_TA
using Test
using Printf
#stucture  contatining molecule's properties: choose molecule from available data
@testset "Molecule properties" begin
    file_name1="data/Molecules/CO2.xlsx"
    @time Atoms,Mol=mol_prop(file_name1)
    ##extraction of properties
    #Number of basis functions Nb
    Nb=Mol.Nb
    #Number of atoms
    Nb_atoms=Mol.Number_atoms
    #Molecule name
    Nb_name=Mol.molecule_name
    #Basis_set
    Basis= Mol.basis_set_name

    #=Atoms is a list of length the number of atoms in the molecule where each atom is a structure
    containing the 3D coordinates in R , exponents and coefficients used for the primitives gaussians=#
    #Extract geometry of atom 1
    Geo=Atoms[1].Geo
    #Extract an orbital from first atoms
    Orbit=Atoms[1].Orbits[1]
    #Each orbital is a structure that contains coeff, expo and name """
    coeff=Orbit.coeff
    expo=Orbit.expo
    name=Orbit.name

    @printf("Molecule properties: \n")
    @printf("Name: %s  Number of atoms: %d  Basis set: %s ", Nb_name, Nb_atoms, Basis)
end


using Test
using NQCBase
using PythonCall

@testset "Single frame" begin
    atoms = Atoms([:H, :C, :O, :N])
    cell = PeriodicCell([10 0 0; 0 10 0; 0 0 10])
    R = rand(3, 4) .* 10

    ase_atoms = convert_to_ase_atoms(atoms, R, cell)
    new_atoms, new_R, new_cell = convert_from_ase_atoms(ase_atoms)
    
    @test new_atoms == atoms
    @test new_cell.vectors ≈ cell.vectors
    @test new_cell.inverse ≈ cell.inverse
    @test new_cell.periodicity ≈ cell.periodicity
    @test new_R ≈ R
end

@testset "Multiple frames" begin
    atoms = Atoms([:H, :C, :O, :N])
    cell = PeriodicCell([10 0 0; 0 10 0; 0 0 10])
    R = [rand(3, 4) .* 10 for _=1:100]

    ase_atoms = convert_to_ase_atoms(atoms, R, cell)
    out = convert_from_ase_atoms.(ase_atoms)
    @test out isa Vector{<:Tuple}
end


using Test
using NQCBase
using PythonCall

@testset "Single frame" begin
    atoms = Atoms([:H, :C, :O, :N])
    cell = PeriodicCell([10 0 0; 0 10 0; 0 0 10])
    R = rand(3, 4) .* 10

    ase_atoms = convert_to_ase_atoms(atoms, R, cell)
    nqcd_structure = convert_from_ase_atoms(ase_atoms)
    
    @test nqcd_structure.atoms == atoms
    @test nqcd_structure.cell.vectors ≈ cell.vectors
    @test nqcd_structure.cell.inverse ≈ cell.inverse
    @test nqcd_structure.cell.periodicity ≈ cell.periodicity
    @test nqcd_structure.positions ≈ R
end

@testset "Multiple frames" begin
    atoms = Atoms([:H, :C, :O, :N])
    cell = PeriodicCell([10 0 0; 0 10 0; 0 0 10])
    R = [rand(3, 4) .* 10 for _=1:100]

    ase_atoms = convert_to_ase_atoms(atoms, R, cell)
    out = convert_from_ase_atoms.(ase_atoms)
    # Correct types
    @test out isa Vector{<:NQCBase.Structure}
    for (frame, positions) in zip(out, R)
        nqcd_structure = frame
        @test nqcd_structure.atoms == atoms
        @test nqcd_structure.cell.vectors ≈ cell.vectors
        @test nqcd_structure.cell.inverse ≈ cell.inverse
        @test nqcd_structure.cell.periodicity ≈ cell.periodicity
        @test nqcd_structure.positions ≈ positions
    end
end

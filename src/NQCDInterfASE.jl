module NQCDInterfASE

using PythonCall
using NQCBase: Atoms, AbstractCell, PeriodicCell, InfiniteCell
import NQCModels
using Unitful, UnitfulAtomic

## IO Interface for ASE

export convert_from_ase_atoms
export convert_to_ase_atoms

const ase = PythonCall.pyimport("ase")

convert_to_ase_atoms(atoms::Atoms, R::Matrix) =
	ase.Atoms(positions=ustrip.(u"Å", R'u"bohr"), symbols=string.(atoms.types))

convert_to_ase_atoms(atoms::Atoms, R::Matrix, ::InfiniteCell) =
	convert_to_ase_atoms(atoms, R)

function convert_to_ase_atoms(atoms::Atoms, R::Matrix, cell::PeriodicCell)
	ase.Atoms(
		positions=ustrip.(u"Å", R'u"bohr"),
		cell=ustrip.(u"Å", cell.vectors'u"bohr"),
		symbols=string.(atoms.types),
		pbc=cell.periodicity)
end

function convert_to_ase_atoms(atoms::Atoms, R::Vector{<:Matrix}, cell::AbstractCell)
	convert_to_ase_atoms.(Ref(atoms), R, Ref(cell))
end

convert_from_ase_atoms(ase_atoms::PythonCall.Py) =
	Atoms(ase_atoms), positions(ase_atoms), Cell(ase_atoms)

Atoms(ase_atoms::PythonCall.Py) = Atoms{Float64}(Symbol.(PythonCall.PyList(ase_atoms.get_chemical_symbols())))

positions(ase_atoms::PythonCall.Py) = austrip.(PythonCall.PyArray(ase_atoms.get_positions())'u"Å")

function Cell(ase_atoms::PythonCall.Py)
	if all(PythonCall.PyArray(ase_atoms.cell.array) .== 0)
		return InfiniteCell()
	else
		return PeriodicCell{Float64}(austrip.(PythonCall.PyArray(ase_atoms.cell.array)'u"Å"), [Bool(x) for x in ase_atoms.pbc])
	end
end

## NQCModels with ASE connection

"""
	AdiabaticASEModel{A} <: AdiabaticModel

Wrapper for an `ase.Atoms` object that has a calculator attached.
This Model will synchronise the positions with the `ase` object and handle the unit conversions.

Implements both `potential` and `derivative!`.

"""
struct ClassicalASEModel{A} <: NQCModels.ClassicalModels.ClassicalModel
	atoms::A
	ndofs::Int
end

export ClassicalASEModel

ClassicalASEModel(atoms::PythonCall.Py) = ClassicalASEModel(atoms, size(pyconvert(Matrix{Float64}, atoms.get_positions()), 2)) # Easy constructor

NQCModels.ndofs(model::ClassicalASEModel) = model.ndofs

function NQCModels.potential(model::ClassicalASEModel, R::AbstractMatrix)
	set_coordinates!(model, R)
	V = model.atoms.get_potential_energy()
	return austrip(pyconvert(eltype(R), V) * u"eV")
end

function NQCModels.derivative!(model::ClassicalASEModel, D::AbstractMatrix, R::AbstractMatrix)
	set_coordinates!(model, R)
	D .= -pyconvert(Matrix{eltype(D)}, model.atoms.get_forces())'
	@. D = austrip(D * u"eV/Å")
	return D
end

function set_coordinates!(model::ClassicalASEModel, R)
	model.atoms.set_positions(ustrip.(auconvert.(u"Å", R')))
end

"""
This module contains methods related to the NQCModels ASE interface that need access to Python types. (e.g. constraint checking)
"""

function NQCModels.mobileatoms(model::ClassicalASEModel, n::Int)
	return symdiff(1:length(model.atoms), [pyconvert(Vector, constraint.get_indices()) .+ 1 for constraint in model.atoms.constraints]...)
end

"""
	ASEFrictionProvider{A} <: TensorialFriction

Obtain the electronic friction from an ASE calculator that implements `get_friction_tensor`.
Assumes that the units of friction are "eV/Å/Å".
Construct by passing the ase atoms object with the calculator already attached.
"""
struct ASEFrictionProvider{A} <: NQCModels.FrictionModels.TensorialFriction
	atoms::A
end

export ASEFrictionProvider

NQCModels.ndofs(::ASEFrictionProvider) = 3

function NQCModels.FrictionModels.get_friction_matrix(model::ASEFrictionProvider, R::AbstractMatrix)
	set_coordinates!(model, R)
	F = pyconvert(Matrix{eltype(R)}, model.atoms.get_friction_tensor()) # Not transposing since the EFT must be symmetric
	@. F = austrip(F * u"eV/Å/Å")
	return F
end

end

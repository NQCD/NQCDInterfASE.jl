module NQCDInterfASE

using PythonCall
using NQCBase: Atoms, AbstractCell, PeriodicCell, Infinitecell
using NQCModels: FrictionModels

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

end

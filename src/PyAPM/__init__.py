"""PyAPM: Python Automated Polymer Modeler for LAMMPS

This package provides automated construction of polymer models (all-atom and coarse-grained)
from SMILES strings, including homopolymers, copolymers, branched, cross-linked, and CG models.
"""

__version__ = "1.0.0"

from .PyAPM import (
	build_oligomer,
	build_homopolymer,
	build_random_copolymer,
	build_alternate_copolymer,
	build_block_copolymer,
	build_sequential_copolymer,
	build_branch_polymer,
	build_crosslink_polymer,
	build_cg_block_copolymer,
)

from .MoleculeBuilder import MoleculeBuilder
from .ChainBuilder import ChainBuilder
from .AtmosphereBuilder import AtmosphereBuilder
from .CrosslinkBuilder import CrosslinkBuilder
from .CGBuilder import Chain, Sequential_chain, Branched_chain

__all__ = [
	"build_oligomer",
	"build_homopolymer",
	"build_random_copolymer",
	"build_alternate_copolymer",
	"build_block_copolymer",
	"build_sequential_copolymer",
	"build_branch_polymer",
	"build_crosslink_polymer",
	"build_cg_block_copolymer",
	"MoleculeBuilder",
	"ChainBuilder",
	"AtmosphereBuilder",
	"CrosslinkBuilder",
	"Chain",
	"Sequential_chain",
	"Branched_chain",
	"__version__",
]
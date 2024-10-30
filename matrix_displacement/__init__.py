"""
Matrix Displacement Method package for structural analysis.
"""

from .core.structure import StructureIdentification
from .core.stiffness import StiffnessComponent, StiffnessMatrixAssembler
from .utils.matrix_utils import build_transformation_matrix, compute_local_stiffness, solve_system

__version__ = '0.1.0'

__all__ = [
    'StructureIdentification',
    'StiffnessComponent',
    'StiffnessMatrixAssembler',
    'build_transformation_matrix',
    'compute_local_stiffness',
    'solve_system'
]

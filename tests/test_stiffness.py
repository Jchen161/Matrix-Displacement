import pytest
import numpy as np
from matrix_displacement.core.stiffness import StiffnessMatrixAssembler, StiffnessComponent
from matrix_displacement.core.structure import StructureIdentification

@pytest.fixture
def simple_structure():
    struct = StructureIdentification()
    n1 = struct.add_node(0, 0)  # Fixed
    n2 = struct.add_node(0, 1, ['ux', 'uy'])  # Free
    e1 = struct.add_element(n1, n2)
    struct.add_force(n2, 'fy', -1000)
    return struct

def test_stiffness_component():
    s1 = StiffnessComponent("EA/L")
    s2 = StiffnessComponent("12EI/L³")
    s3 = s1 + s2
    assert str(s3) == "EA/L + 12EI/L³"
    
    s4 = StiffnessComponent("0")
    s5 = s4 + s1
    assert str(s5) == "EA/L"

def test_assembler_initialization(simple_structure):
    assembler = StiffnessMatrixAssembler(simple_structure)
    assert assembler.structure == simple_structure
    assert len(assembler.element_matrices) == 0

def test_element_properties(simple_structure):
    assembler = StiffnessMatrixAssembler(simple_structure)
    assembler.set_default_properties(E=200e9, A=0.001, I=1e-6)
    props = assembler.get_element_properties(0)
    assert props['E'] == 200e9
    assert props['A'] == 0.001
    assert props['I'] == 1e-6

def test_element_dofs(simple_structure):
    assembler = StiffnessMatrixAssembler(simple_structure)
    active_dofs, indices = assembler.get_element_dofs(0)
    assert len(active_dofs) == 2  # ux and uy for free node
    assert all(idx is not None for idx in indices)

def test_stiffness_computation(simple_structure):
    assembler = StiffnessMatrixAssembler(simple_structure)
    assembler.set_default_properties(E=200e9, A=0.001, I=1e-6)
    k_local_sym, k_local_num, k_global_sym, k_global_num = assembler.compute_element_stiffness(0)
    
    assert k_local_sym is not None
    assert k_local_num is not None
    assert k_global_sym is not None
    assert k_global_num is not None
    
    # Test matrix dimensions
    n_dofs = len(assembler.element_matrices[0]['active_dofs'])
    assert k_local_num.shape == (n_dofs, n_dofs)
    assert k_global_num.shape == (n_dofs, n_dofs)

def test_global_assembly(simple_structure):
    assembler = StiffnessMatrixAssembler(simple_structure)
    assembler.set_default_properties(E=200e9, A=0.001, I=1e-6)
    assembler.compute_element_stiffness(0)
    K_sym, K_num = assembler.assemble_global_stiffness_matrix()
    
    assert K_sym is not None
    assert K_num is not None
    assert K_num.shape == (simple_structure.dof_count, simple_structure.dof_count)

def test_displacement_solution(simple_structure):
    assembler = StiffnessMatrixAssembler(simple_structure)
    assembler.set_default_properties(E=200e9, A=0.001, I=1e-6)
    assembler.compute_element_stiffness(0)
    assembler.assemble_global_stiffness_matrix()
    d_sym, d_num = assembler.solve_displacements()
    
    assert d_sym is not None
    assert d_num is not None
    assert len(d_num) == simple_structure.dof_count

if __name__ == '__main__':
    pytest.main([__file__])

import pytest
import numpy as np
from matrix_displacement.core.structure import StructureIdentification

def test_node_creation():
    struct = StructureIdentification()
    
    # Test basic node creation
    node_id = struct.add_node(0, 1)
    assert node_id == 0
    assert np.array_equal(struct.nodes[0]['coordinates'], np.array([0, 1]))
    assert struct.node_count == 1
    
    # Test node with DOFs
    node_id = struct.add_node(1, 1, ['ux', 'uy'])
    assert node_id == 1
    assert struct.nodes[1]['has_dof']['ux'] == True
    assert struct.nodes[1]['has_dof']['uy'] == True
    assert struct.nodes[1]['has_dof']['theta'] == False

def test_element_creation():
    struct = StructureIdentification()
    n1 = struct.add_node(0, 0)
    n2 = struct.add_node(1, 0)
    
    # Test element creation
    e1 = struct.add_element(n1, n2)
    assert e1 == 0
    assert struct.elements[0]['nodes'] == (n1, n2)
    assert struct.elements[0]['length'] == pytest.approx(1.0)
    assert struct.elements[0]['cos_theta'] == pytest.approx(1.0)
    assert struct.elements[0]['sin_theta'] == pytest.approx(0.0)

def test_force_application():
    struct = StructureIdentification()
    n1 = struct.add_node(0, 0, ['ux', 'uy'])
    
    # Test force application
    struct.add_force(n1, 'fx', 1000)
    assert struct.forces[struct.nodes[n1]['dofs']['ux']] == 1000
    
    # Test invalid force application
    with pytest.raises(ValueError):
        struct.add_force(n1, 'M', 100)  # No rotational DOF

def test_shared_dofs():
    struct = StructureIdentification()
    n1 = struct.add_node(0, 0, ['ux', 'uy'])
    n2 = struct.add_node(1, 0, dofs=None, 
                        shared_dofs={'ux': (n1, 'ux')})
    
    # Test shared DOF
    assert struct.nodes[n2]['dofs']['ux'] == struct.nodes[n1]['dofs']['ux']
    assert struct.nodes[n2]['has_dof']['ux'] == True

def test_force_vector():
    struct = StructureIdentification()
    n1 = struct.add_node(0, 0, ['ux', 'uy'])
    struct.add_force(n1, 'fx', 1000)
    
    F_sym, F_num = struct.get_force_vector()
    assert len(F_sym) == struct.dof_count
    assert len(F_num) == struct.dof_count
    assert F_num[0] == 1000  # First DOF should have the force

if __name__ == '__main__':
    pytest.main([__file__])

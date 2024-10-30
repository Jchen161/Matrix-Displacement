import pytest
import numpy as np
from matrix_displacement.utils.matrix_utils import (
    build_transformation_matrix,
    compute_local_stiffness,
    solve_system
)

def test_transformation_matrix():
    active_dofs = [(0, 'ux'), (0, 'uy')]
    cos_theta = 1.0  # 0 degrees
    sin_theta = 0.0
    
    T = build_transformation_matrix(active_dofs, cos_theta, sin_theta)
    assert T.shape == (2, 2)
    assert np.allclose(T, np.eye(2))

def test_compute_local_stiffness():
    # Test parameters
    node_i, node_j = 0, 1
    E = 200e9  # Steel
    A = 0.001  # 1000 mm²
    I = 1e-6   # Moment of inertia
    L = 1.0    # 1m length
    
    # Test axial terms
    k_axial = compute_local_stiffness(node_i, node_j, 'ux', 'ux', E, A, I, L)
    assert k_axial == -E * A / L  # Different nodes, negative value
    
    # Test lateral terms
    k_lateral = compute_local_stiffness(node_i, node_i, 'uy', 'uy', E, A, I, L)
    assert k_lateral == 12 * E * I / (L**3)  # Same node, positive value
    
    # Test rotational terms
    k_rot = compute_local_stiffness(node_i, node_j, 'theta', 'theta', E, A, I, L)
    assert k_rot == 2 * E * I / L  # Different nodes
    
    # Test coupling terms
    k_coupling = compute_local_stiffness(node_i, node_j, 'uy', 'theta', E, A, I, L)
    assert k_coupling == 6 * E * I / (L**2)  # Coupling term

def test_invalid_properties():
    node_i, node_j = 0, 1
    result = compute_local_stiffness(node_i, node_j, 'ux', 'ux', None, 0.001, 1e-6, 1.0)
    assert result is None

def test_solve_system():
    K = np.array([[2, -1], [-1, 2]])
    F = np.array([1, 0])
    
    d = solve_system(K, F)
    assert len(d) == 2
    assert np.allclose(K @ d, F)

def test_singular_matrix():
    K = np.array([[1, 1], [1, 1]])  # Singular matrix
    F = np.array([1, 1])
    
    result = solve_system(K, F)
    assert result is None

if __name__ == '__main__':
    pytest.main([__file__])
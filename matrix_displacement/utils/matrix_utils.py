"""
Utility functions for matrix operations in structural analysis.
Provides core mathematical operations for the matrix displacement method.
"""

import numpy as np
from typing import List, Tuple, Optional, Union

def build_transformation_matrix(active_dofs: List[Tuple[int, str]], 
                              cos_theta: float, 
                              sin_theta: float) -> np.ndarray:
    """
    Build coordinate transformation matrix for element local to global conversion.

    Args:
        active_dofs: List of active degrees of freedom [(node_number, dof_type), ...]
        cos_theta: Cosine of element angle with global x-axis
        sin_theta: Sine of element angle with global x-axis

    Returns:
        np.ndarray: Transformation matrix
    """
    # Complete transformation matrix template
    T_general = np.array([
        [cos_theta,  sin_theta, 0, 0, 0, 0],
        [-sin_theta, cos_theta, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, cos_theta,  sin_theta, 0],
        [0, 0, 0, -sin_theta, cos_theta, 0],
        [0, 0, 0, 0, 0, 1]
    ])
    
    n_dofs = len(active_dofs)
    T = np.zeros((n_dofs, n_dofs))
    
    # Map DOF types to indices in the transformation matrix
    dof_map = {
        ('ux', 0): 0,    # First node x direction
        ('uy', 0): 1,    # First node y direction
        ('theta', 0): 2, # First node rotation
        ('ux', 1): 3,    # Second node x direction
        ('uy', 1): 4,    # Second node y direction
        ('theta', 1): 5  # Second node rotation
    }
    
    # Get unique nodes and create mapping
    node_nums = sorted(list(set(node for node, _ in active_dofs)))
    node_map = {node: idx for idx, node in enumerate(node_nums)}
    
    # Build transformation matrix based on active DOFs
    for i, (node_i, dof_i) in enumerate(active_dofs):
        for j, (node_j, dof_j) in enumerate(active_dofs):
            row = dof_map[(dof_i, node_map[node_i])]
            col = dof_map[(dof_j, node_map[node_j])]
            T[i, j] = T_general[row, col]
    
    return T

def compute_local_stiffness(node_i: int, node_j: int, dof_i: str, dof_j: str, 
                          E: float, A: float, I: float, L: float) -> float:
    """
    Compute a single term of the element stiffness matrix in local coordinates.
    
    Args:
        node_i: First node number
        node_j: Second node number
        dof_i: First degree of freedom type ('ux', 'uy', 'theta')
        dof_j: Second degree of freedom type ('ux', 'uy', 'theta')
        E: Elastic modulus
        A: Cross-sectional area
        I: Moment of inertia
        L: Element length
    
    Returns:
        float: Stiffness matrix term
    """
    if not all(x is not None for x in [E, A, I]):
        return None
    
    same_node = node_i == node_j
    
    # Precalculate common values
    EA_L = E * A / L
    EI_L = E * I / L
    EI_L2 = EI_L / L
    EI_L3 = EI_L2 / L
    
    # Axial terms (ux-ux)
    if dof_i == 'ux' and dof_j == 'ux':
        return EA_L if same_node else -EA_L
    
    # Lateral terms (uy-uy)
    elif dof_i == 'uy' and dof_j == 'uy':
        return 12 * EI_L3 if same_node else -12 * EI_L3
    
    # Rotational terms (theta-theta)
    elif dof_i == 'theta' and dof_j == 'theta':
        return 4 * EI_L if same_node else 2 * EI_L
    
    # Lateral displacement-rotation coupling terms
    elif (dof_i == 'uy' and dof_j == 'theta') or (dof_i == 'theta' and dof_j == 'uy'):
        if same_node:
            return 6 * EI_L2 if dof_i == 'uy' else -6 * EI_L2
        else:
            if (dof_i == 'uy' and node_i < node_j) or (dof_j == 'uy' and node_j < node_i):
                return 6 * EI_L2
            else:
                return -6 * EI_L2
    
    # All other combinations are zero
    return 0.0

def solve_system(K: np.ndarray, F: np.ndarray) -> Optional[np.ndarray]:
    """
    Solve the system of equations K∆ = F.

    Args:
        K: Global stiffness matrix
        F: Force vector

    Returns:
        Optional[np.ndarray]: Displacement vector or None if system is singular

    Note:
        Uses numpy's linear algebra solver with error handling for singular matrices
    """
    try:
        # Check matrix condition
        if np.linalg.cond(K) > 1e15:
            print("Warning: System is ill-conditioned")
            
        # Solve system
        displacements = np.linalg.solve(K, F)
        return displacements
        
    except np.linalg.LinAlgError as e:
        print(f"Error solving system: {e}")
        print("Check if structure is properly constrained")
        return None

def calculate_element_forces(k_local: np.ndarray, 
                           T: np.ndarray, 
                           d_global: np.ndarray) -> np.ndarray:
    """
    Calculate element end forces in local coordinates.

    Args:
        k_local: Element stiffness matrix in local coordinates
        T: Transformation matrix
        d_global: Global displacement vector

    Returns:
        np.ndarray: Element end forces in local coordinates
    """
    # Transform global displacements to local
    d_local = T @ d_global
    
    # Calculate forces
    f_local = k_local @ d_local
    
    return f_local
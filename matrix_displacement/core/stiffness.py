"""
Stiffness matrix computation and assembly module.
Provides classes and methods for calculating and managing structural stiffness matrices.
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Union
from ..utils.matrix_utils import (
    build_transformation_matrix,
    compute_local_stiffness,
    solve_system
)

class StiffnessComponent:
    """Handles symbolic expressions for stiffness matrix elements."""
    
    def __init__(self, symbolic: str):
        self.symbolic = symbolic
    
    def __repr__(self) -> str:
        return f"{self.symbolic}"
    
    def __add__(self, other: 'StiffnessComponent') -> 'StiffnessComponent':
        if isinstance(other, StiffnessComponent):
            if self.symbolic == "0":
                return other
            if other.symbolic == "0":
                return self
            return StiffnessComponent(f"{self.symbolic} + {other.symbolic}")
        return self

    def __radd__(self, other: 'StiffnessComponent') -> 'StiffnessComponent':
        return self.__add__(other)

class StiffnessMatrixAssembler:
    """Assembles and manages structural stiffness matrices."""
    
    def __init__(self, structure):
        self.structure = structure
        self.element_matrices: Dict = {}
        self.element_properties: Dict = {}
        self.global_matrix_symbolic = None
        self.global_matrix_numeric = None
        self.displacements_symbolic = None
        self.displacements_numeric = None
        self.default_properties = {
            'E': None,
            'A': None,
            'I': None
        }
    
    def set_element_property(self, element_id: int, E: Optional[float] = None, 
                           A: Optional[float] = None, I: Optional[float] = None):
        """Set properties for a specific element."""
        if element_id not in self.structure.elements:
            raise ValueError(f"Element {element_id} does not exist")
        self.element_properties[element_id] = {'E': E, 'A': A, 'I': I}
    
    def set_default_properties(self, E: Optional[float] = None, 
                             A: Optional[float] = None, I: Optional[float] = None):
        """Set default properties for all elements."""
        self.default_properties = {'E': E, 'A': A, 'I': I}
    
    def get_element_properties(self, element_id: int) -> Dict[str, float]:
        """Get properties for an element, using defaults if not specifically set."""
        properties = {}
        element_props = self.element_properties.get(element_id, {})
        
        for prop in ['E', 'A', 'I']:
            element_value = element_props.get(prop)
            properties[prop] = (element_value if element_value is not None 
                              else self.default_properties[prop])
        return properties
    
    def get_element_dofs(self, element_id: int) -> Tuple[List[Tuple[int, str]], List[int]]:
        """Get element's degrees of freedom information."""
        element = self.structure.elements[element_id]
        node_i, node_j = element['nodes']
        
        active_dofs = []
        dof_indices = []
        
        for node_id in [node_i, node_j]:
            node = self.structure.nodes[node_id]
            for dof_type in ['ux', 'uy', 'theta']:
                if node['has_dof'][dof_type]:
                    active_dofs.append((node_id, dof_type))
                    dof_indices.append(node['dofs'][dof_type])
        
        return active_dofs, dof_indices
    
    def compute_element_stiffness(self, element_id: int) -> Tuple:
        """Compute element stiffness matrices."""
        properties = self.get_element_properties(element_id)
        E, A, I = properties['E'], properties['A'], properties['I']
        element = self.structure.elements[element_id]
        L = element['length']
        c = element['cos_theta']
        s = element['sin_theta']
        
        active_dofs, dof_indices = self.get_element_dofs(element_id)
        n_dofs = len(active_dofs)
        
        # Compute symbolic matrices
        k_local_symbolic = np.empty((n_dofs, n_dofs), dtype=object)
        k_local_numeric = np.zeros((n_dofs, n_dofs)) if all(x is not None for x in [E, A, I]) else None
        
        for i, (node_i, dof_i) in enumerate(active_dofs):
            for j, (node_j, dof_j) in enumerate(active_dofs):
                k_local_symbolic[i, j] = self.get_symbolic_stiffness_term(
                    node_i, dof_i, node_j, dof_j, element_id
                )
                if k_local_numeric is not None:
                    k_local_numeric[i, j] = compute_local_stiffness(
                        node_i, node_j, dof_i, dof_j, E, A, I, L
                    )
        
        # Create transformation matrix
        T = build_transformation_matrix(active_dofs, c, s)
        
        # Compute global matrices
        k_global_symbolic = np.zeros((n_dofs, n_dofs), dtype=object)
        k_global_numeric = None
        
        # Symbolic global matrix
        for i in range(n_dofs):
            for j in range(n_dofs):
                k_global_symbolic[i, j] = StiffnessComponent("0")
                for k in range(n_dofs):
                    for l in range(n_dofs):
                        if T[k,i] != 0 and T[l,j] != 0 and str(k_local_symbolic[k,l]) != "0":
                            term = self._create_symbolic_term(k_local_symbolic[k,l], T[k,i], T[l,j])
                            k_global_symbolic[i, j] = k_global_symbolic[i, j] + term
        
        # Numeric global matrix
        if k_local_numeric is not None:
            k_global_numeric = T.T @ k_local_numeric @ T
        
        # Store matrices
        self.element_matrices[element_id] = {
            'active_dofs': active_dofs,
            'dof_indices': dof_indices,
            'local_symbolic': k_local_symbolic,
            'local_numeric': k_local_numeric,
            'global_symbolic': k_global_symbolic,
            'global_numeric': k_global_numeric,
            'transformation': T,
            'parameters': {'E': E, 'A': A, 'I': I, 'L': L} if k_local_numeric is not None else None
        }
        
        return k_local_symbolic, k_local_numeric, k_global_symbolic, k_global_numeric

    def _create_symbolic_term(self, base_term: StiffnessComponent, 
                            t1: float, t2: float) -> StiffnessComponent:
        """Create symbolic term with transformation coefficients."""
        if t1 == 1 and t2 == 1:
            return base_term
        elif t1 == 1:
            return StiffnessComponent(f"{base_term}×{t2}")
        elif t2 == 1:
            return StiffnessComponent(f"{t1}×{base_term}")
        return StiffnessComponent(f"{t1}×{base_term}×{t2}")

    def get_symbolic_stiffness_term(self, node_i: int, dof_i: str, 
                                  node_j: int, dof_j: str, 
                                  element_id: int) -> StiffnessComponent:
        """Get symbolic term for stiffness matrix."""
        same_node = node_i == node_j
        e_suffix = f"_{element_id}"
        
        if dof_i == 'ux' and dof_j == 'ux':
            return StiffnessComponent(f"EA{e_suffix}/l{e_suffix}") if same_node else \
                   StiffnessComponent(f"-EA{e_suffix}/l{e_suffix}")
        
        elif dof_i == 'uy' and dof_j == 'uy':
            return StiffnessComponent(f"12EI{e_suffix}/l{e_suffix}³") if same_node else \
                   StiffnessComponent(f"-12EI{e_suffix}/l{e_suffix}³")
        
        elif dof_i == 'theta' and dof_j == 'theta':
            return StiffnessComponent(f"4EI{e_suffix}/l{e_suffix}") if same_node else \
                   StiffnessComponent(f"2EI{e_suffix}/l{e_suffix}")
        
        elif (dof_i == 'uy' and dof_j == 'theta') or (dof_i == 'theta' and dof_j == 'uy'):
            if same_node:
                return StiffnessComponent(f"6EI{e_suffix}/l{e_suffix}²") if dof_i == 'uy' else \
                       StiffnessComponent(f"-6EI{e_suffix}/l{e_suffix}²")
            else:
                if (dof_i == 'uy' and node_i < node_j) or (dof_j == 'uy' and node_j < node_i):
                    return StiffnessComponent(f"6EI{e_suffix}/l{e_suffix}²")
                else:
                    return StiffnessComponent(f"-6EI{e_suffix}/l{e_suffix}²")
        
        return StiffnessComponent("0")

    def assemble_global_stiffness_matrix(self) -> Tuple:
        """Assemble global stiffness matrix."""
        n_dofs = self.structure.dof_count
        K_sym = [[StiffnessComponent("0") for _ in range(n_dofs)] 
                for _ in range(n_dofs)]
        
        has_numeric = all(matrices.get('global_numeric') is not None 
                         for matrices in self.element_matrices.values())
        K_num = np.zeros((n_dofs, n_dofs)) if has_numeric else None
        
        for element_id, matrices in self.element_matrices.items():
            k_global_sym = matrices['global_symbolic']
            dof_indices = matrices['dof_indices']
            
            for i, gi in enumerate(dof_indices):
                for j, gj in enumerate(dof_indices):
                    if gi is not None and gj is not None:
                        K_sym[gi][gj] = K_sym[gi][gj] + k_global_sym[i][j]
            
            if has_numeric:
                k_global_num = matrices['global_numeric']
                for i, gi in enumerate(dof_indices):
                    for j, gj in enumerate(dof_indices):
                        if gi is not None and gj is not None:
                            K_num[gi, gj] += k_global_num[i, j]
        
        self.global_matrix_symbolic = K_sym
        self.global_matrix_numeric = K_num
        return K_sym, K_num

    def solve_displacements(self) -> Tuple:
        """Solve for displacements using K∆ = F."""
        if self.global_matrix_symbolic is None:
            raise ValueError("Global stiffness matrix not yet assembled")
        
        F_sym, F_num = self.structure.get_force_vector()
        n_dofs = len(self.global_matrix_symbolic)
        
        self.displacements_symbolic = [
            StiffnessComponent(f"∆_{i}") for i in range(n_dofs)
        ]
        
        if (self.global_matrix_numeric is not None and 
            len(self.structure.forces) > 0):
            try:
                self.displacements_numeric = solve_system(
                    self.global_matrix_numeric, F_num
                )
                return self.displacements_symbolic, self.displacements_numeric
            except Exception as e:
                print(f"Error solving numeric equation: {e}")
                return self.displacements_symbolic, None
        
        return self.displacements_symbolic, None

    def print_solution(self):
        """Print the solution in a readable format."""
        if self.displacements_symbolic is None:
            print("No solution available. Please solve the displacement equations first.")
            return
        
        print("\nDisplacement Solution:")
        print("-" * 50)
        
        node_displacements = {}
        for node_id, node in self.structure.nodes.items():
            node_displacements[node_id] = {
                'ux': None, 'uy': None, 'theta': None,
                'ux_num': None, 'uy_num': None, 'theta_num': None
            }
            
            for dof_type in ['ux', 'uy', 'theta']:
                dof_num = node['dofs'][dof_type]
                if dof_num is not None:
                    node_displacements[node_id][dof_type] = \
                        self.displacements_symbolic[dof_num]
                    if self.displacements_numeric is not None:
                        node_displacements[node_id][f'{dof_type}_num'] = \
                            self.displacements_numeric[dof_num]
        
        for node_id, displacements in node_displacements.items():
            print(f"\nNode {node_id}:")
            for dof_type in ['ux', 'uy', 'theta']:
                symbol = displacements[dof_type]
                numeric = displacements[f'{dof_type}_num']
                
                if symbol is not None:
                    print(f"  {dof_type}: {symbol}", end="")
                    if numeric is not None:
                        print(f" = {numeric:.6e}")
                    else:
                        print(" (no numeric solution)")
                else:
                    print(f"  {dof_type}: Constrained")
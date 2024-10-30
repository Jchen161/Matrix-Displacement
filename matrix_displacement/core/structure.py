"""
Structure identification and management module.
Handles node and element definitions, degrees of freedom, and force applications.
"""

import numpy as np
from typing import List, Dict, Optional, Tuple, Union

class StructureIdentification:
    def __init__(self):
        """Initialize a new structure with empty nodes, elements, and forces."""
        self.nodes: Dict = {}
        self.elements: Dict = {}
        self.node_count: int = 0
        self.element_count: int = 0
        self.forces: Dict = {}
        self.dof_count: int = 0
    
    def add_node(self, x: float, y: float, 
                 dofs: Optional[List[str]] = None,
                 shared_dofs: Optional[Dict] = None) -> int:
        """
        Add a node with specified degrees of freedom and shared DOFs.

        Args:
            x: X-coordinate of the node
            y: Y-coordinate of the node
            dofs: List of degrees of freedom ['ux', 'uy', 'theta']
            shared_dofs: Dictionary specifying shared DOFs with other nodes

        Returns:
            int: Node identifier
        """
        node_id = self.node_count
        
        node_dofs = {
            'ux': None,      # x-direction displacement
            'uy': None,      # y-direction displacement
            'theta': None    # rotation
        }
        
        has_dof = {
            'ux': False,
            'uy': False,
            'theta': False
        }
        
        if shared_dofs:
            for dof_type, (other_node, other_dof) in shared_dofs.items():
                if other_node in self.nodes and dof_type in node_dofs:
                    shared_dof_num = self.nodes[other_node]['dofs'][other_dof]
                    if shared_dof_num is not None:
                        node_dofs[dof_type] = shared_dof_num
                        has_dof[dof_type] = True
        
        if dofs:
            for dof_type in dofs:
                if dof_type in node_dofs and not has_dof[dof_type]:
                    node_dofs[dof_type] = self.dof_count
                    has_dof[dof_type] = True
                    self.dof_count += 1
        
        self.nodes[node_id] = {
            'coordinates': np.array([x, y]),
            'dofs': node_dofs,
            'has_dof': has_dof
        }
        
        self.node_count += 1
        return node_id

    def add_element(self, node_i: int, node_j: int) -> int:
        """
        Add an element connecting two nodes.

        Args:
            node_i: First node ID
            node_j: Second node ID

        Returns:
            int: Element identifier
        
        Raises:
            ValueError: If either node does not exist
        """
        if node_i not in self.nodes or node_j not in self.nodes:
            raise ValueError("Both nodes must exist before creating an element")
            
        element_id = self.element_count
        
        coord_i = self.nodes[node_i]['coordinates']
        coord_j = self.nodes[node_j]['coordinates']
        
        dx = coord_j[0] - coord_i[0]
        dy = coord_j[1] - coord_i[1]
        length = np.sqrt(dx**2 + dy**2)
        cos_theta = dx/length
        sin_theta = dy/length
        
        element_dofs = {
            'start_node': {
                'ux': self.nodes[node_i]['dofs']['ux'],
                'uy': self.nodes[node_i]['dofs']['uy'],
                'theta': self.nodes[node_i]['dofs']['theta']
            },
            'end_node': {
                'ux': self.nodes[node_j]['dofs']['ux'],
                'uy': self.nodes[node_j]['dofs']['uy'],
                'theta': self.nodes[node_j]['dofs']['theta']
            }
        }
        
        self.elements[element_id] = {
            'nodes': (node_i, node_j),
            'length': length,
            'cos_theta': cos_theta,
            'sin_theta': sin_theta,
            'dofs': element_dofs
        }
        
        self.element_count += 1
        return element_id

    def add_force(self, node_id: int, force_type: str, magnitude: float) -> None:
        """
        Add force or moment to a node.

        Args:
            node_id: Node identifier
            force_type: 'fx' (x-direction force), 'fy' (y-direction force), 
                       or 'M' (moment)
            magnitude: Force/moment magnitude
        
        Raises:
            ValueError: If node doesn't exist or has no corresponding DOF
        """
        if node_id not in self.nodes:
            raise ValueError(f"Node {node_id} does not exist")
            
        force_to_dof = {'fx': 'ux', 'fy': 'uy', 'M': 'theta'}
        dof_type = force_to_dof[force_type]
        
        dof_num = self.nodes[node_id]['dofs'][dof_type]
        
        if dof_num is None:
            raise ValueError(f"Node {node_id} has no {dof_type} degree of freedom")
            
        self.forces[dof_num] = magnitude

    def get_force_vector(self) -> Tuple[List, np.ndarray]:
        """
        Create the global force vector.

        Returns:
            Tuple containing:
            - List of symbolic force components
            - Numpy array of numeric force values
        """
        from .stiffness import StiffnessComponent
        
        F_sym = [StiffnessComponent("0") for _ in range(self.dof_count)]
        F_num = np.zeros(self.dof_count)
        
        for dof_num, force in self.forces.items():
            F_sym[dof_num] = StiffnessComponent(f"F_{dof_num}")
            F_num[dof_num] = force
            
        return F_sym, F_num

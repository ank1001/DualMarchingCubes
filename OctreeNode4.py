import numpy as np
from FeatureIsolation4 import Feature_Isolation, compute_f
from ModMarchingCubes2 import mod_marching_cubes
import mcubes


class OctreeNode:


    def __init__(self, base_coords, side_length, number_of_points, is_leaf=False) -> None:
        self.base_coords = base_coords
        self.side_length = side_length
        self.number_of_points = number_of_points
        self.QEF_threshold = 0.01
        self.is_leaf = is_leaf
        self.cell_params = None #[w, x, y, z]
        self.QEF = None
        self.children = [None]*8
        self.get_QEF()

    def __str__(self):
        return (f"OctreeNode("
            f"base_coords={self.base_coords}, "
            f"side_length={self.side_length}, "
            f"number_of_points={self.number_of_points}, "
            f"cell_params={self.cell_params}, "
            f"QEF={self.QEF}, "
            f"is_leaf={self.is_leaf}, "
            f"children={['present' if child else 'None' for child in self.children]})")

    def get_QEF(self):
        obj = Feature_Isolation(self.base_coords, self.side_length, field_values, isolevel, self.number_of_points)
        self.cell_params=obj.cell_params
        self.QEF = obj.QEF
        if self.QEF<self.QEF_threshold:
            self.is_leaf  = True
    
    def subdivide(self):
        if self.is_leaf:
            return
        half_size = self.side_length / 2
        offsets = [
            [0, 0, 0],
            [half_size, 0, 0],
            [0, half_size, 0],
            [half_size, half_size, 0],
            [0, 0, half_size],
            [half_size, 0, half_size],
            [0, half_size, half_size],
            [half_size, half_size, half_size],
        ]
        new_number_of_points = self.number_of_points//8
        
        is_child_leaf = False
        if new_number_of_points==0:
            is_child_leaf=True
        for i in range(8):
            child_coords = self.base_coords + offsets[i]
            self.children[i] = OctreeNode(child_coords, half_size, new_number_of_points, is_child_leaf)


    def get_or_copy(self):
        if self.is_leaf():
            return [self]*8
        return self.children
    

def faceProc(node):
    # print(f'FACEPROC')
    if node.is_leaf:
        return
    node.subdivide()

    # recursion on its 8 children
    for child in node.children:
        if child:
            faceProc(child)

    # Recursive calls to edgeProc for edge-adjacent pairs
    edge_pairs = [(node.children[0], node.children[1]),
                  (node.children[2], node.children[3]),
                  (node.children[4], node.children[5]),
                  (node.children[6], node.children[7]),
                  
                  (node.children[0], node.children[4]),
                  (node.children[1], node.children[5]),
                  (node.children[2], node.children[6]),
                  (node.children[3], node.children[7]),]
    
    for child1, child2 in edge_pairs:
        edgeProc(child1, child2)

    vertProc(node.children[0], node.children[1], node.children[2], node.children[3], node.children[4], node.children[5], node.children[6], node.children[7])


def edgeProc(node1, node2):
    # print(f'EDGEPROC')
    if not node1 or not node2:
        return
    
    if node1.is_leaf and node2.is_leaf:
        return
    
    if (node1.base_coords[0] != node2.base_coords[0]):
        edge_pairs = [(node1.children[1], node2.children[0]),
                    (node1.children[3], node2.children[2]),
                    (node1.children[5], node2.children[4]),
                    (node1.children[7], node2.children[6])]
        for child1, child2 in edge_pairs:
            edgeProc(child1, child2)
        vertProc(node1.children[1], node2.children[0], node1.children[3], node2.children[2], node1.children[5], node2.children[4], node1.children[7], node2.children[6])
    else:
        edge_pairs = [(node1.children[4], node2.children[0]),
                    (node1.children[5], node2.children[1]),
                    (node1.children[6], node2.children[2]),
                    (node1.children[7], node2.children[3])]
        for child1, child2 in edge_pairs:
            edgeProc(child1, child2)

        vertProc(node1.children[4], node1.children[5], node1.children[6], node1.children[7], node2.children[0], node2.children[1], node2.children[2], node2.children[3])

def vertProc(node1, node2, node3, node4, node5, node6, node7, node8):
    # print(f'VERTPROC')
    if not node1 or not node2 or not node3 or not node4 or not node4 or not node5 or not node6 or not node7 or not node8:
        return

    if all([node1.is_leaf, node2.is_leaf, node3.is_leaf, node4.is_leaf, node5.is_leaf, node6.is_leaf, node7.is_leaf, node8.is_leaf]):
        # get the representative vertex of each cell and for a grid
        cell_params = np.array([
            node1.cell_params,
            node2.cell_params,
            node3.cell_params,
            node4.cell_params,
            node5.cell_params,
            node6.cell_params,
            node7.cell_params,
            node8.cell_params
        ])
        # print(f'cellparams:{cell_params}')
        verts, fcs = mod_marching_cubes(cell_params, isolevel, field_values)
        # print(f'verts:{verts}, fcs:{fcs}')
        dmcoutput.insert(verts, fcs)
        # dual_cells.append(vertices)
        return
    vertProc(node1.children[7], node1.children[6], node1.children[5], node1.children[4], node1.children[3], node1.children[2], node1.children[1], node1.children[0], )

class DMCOutput:
    def __init__(self) -> None:
        self.vertices = []
        self.triangles =[]
        self.vertices_len = 0

    def insert(self, verts, fcs):
        if fcs is not None:
            v = verts.tolist()
            f = fcs.tolist()
            for vert in v:
                self.vertices.append(vert)
            for fc in f:
                fc[0]+= self.vertices_len
                fc[1]+= self.vertices_len
                fc[2]+= self.vertices_len
                self.triangles.append(fc)
            self.vertices_len += len(v)


def dual_marching_cubes(fld_values, isol):
    global field_values
    field_values = fld_values
    global isolevel
    global dmcoutput
    dmcoutput = DMCOutput()
    isolevel = isol
    n = field_values.shape[0]
    n = (n-1)**3
    node = OctreeNode(np.array([0,0,0]), 1, n)
    faceProc(node)
    vertices = np.array(dmcoutput.vertices)
    triangles = np.array(dmcoutput.triangles)
    return vertices, triangles

def main():
    density_threshold = 0.3
    X, Y, Z = np.mgrid[:10, :10, :10]
    u =  (X-4)**2 + (Y-4)**2 + (Z-4)**2 - 3.8**2
    maxu = np.max(u)
    u = u/maxu
    print('Running dual marching cubes ... ')
    vertices, triangles = dual_marching_cubes(u, density_threshold)
    # print(f'vertices:{vertices}, triangles:{triangles}')
    mcubes.export_obj(vertices, triangles, 'mydual_sphere.obj')
    print('exported to obj file.')

main()
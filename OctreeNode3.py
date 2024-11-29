import numpy as np
import math
import random
import mcubes
from FeatureIsolation3 import Feature_Isolation, compute_f


class OctreeNode:


    def __init__(self, base_coords, side_length, parent_QEF = 10) -> None:
        self.base_coords = base_coords
        self.side_length = side_length
        self.isolevel = isolevel
        self.QEF = self.compute_QEF()
        self.QEF_threshold = 0.00001
        self.children = [None]*8
        self.parent_QEF = parent_QEF

    def compute_QEF(self):
        obj = Feature_Isolation(self.base_coords, self.side_length, field_values, self.isolevel)
        return obj.QEF
    
    def is_leaf(self):
        # if abs(self.QEF- self.QEF_threshold)>2:
        if self.QEF> self.parent_QEF:
            return True
        else:
            # create children
            # print(f'create children')
            self.children[0] = OctreeNode(self.base_coords, self.side_length/2.0, self.QEF)
            self.children[1] = OctreeNode(self.base_coords + np.array([self.side_length/2, 0, 0]), self.side_length/2, self.QEF)
            self.children[2] = OctreeNode(self.base_coords + np.array([0, self.side_length/2, 0]), self.side_length/2.0, self.QEF)
            self.children[3] = OctreeNode(self.base_coords + np.array([self.side_length/2, self.side_length/2, 0]), self.side_length/2.0, self.QEF)

            self.children[4] = OctreeNode(self.base_coords + np.array([0, 0, self.side_length/2]), self.side_length/2.0, self.QEF)
            self.children[5] = OctreeNode(self.base_coords + np.array([self.side_length/2, 0, self.side_length/2]), self.side_length/2.0, self.QEF)
            self.children[6] = OctreeNode(self.base_coords + np.array([0, self.side_length/2, self.side_length/2]), self.side_length/2.0, self.QEF)
            self.children[7] = OctreeNode(self.base_coords + np.array([self.side_length/2, self.side_length/2, self.side_length/2]), self.side_length/2.0, self.QEF)

            return False

    def get_or_copy(self):
        if self.is_leaf():
            return [self]*8
        return self.children
    



def faceProc(q0):
    if q0.is_leaf():
        return
    
    # Recrsive calls to faceProc for the children
    children = q0.get_or_copy()
    for child in children:
        # print(f'{child.base_coords}, QEF: {child.QEF}')
        faceProc(child)

    # Recursive calls to edgeProc for edge-adjacent pairs
    edge_pairs = [(children[0], children[1]),
                  (children[2], children[3]),
                  (children[4], children[5]),
                  (children[6], children[7]),
                  
                  (children[0], children[4]),
                  (children[1], children[5]),
                  (children[2], children[6]),
                  (children[3], children[7]),]
    
    for child1, child2 in edge_pairs:
        edgeProc(child1, child2)

    vertProc(children[0], children[1], children[2], children[3], children[4], children[5], children[6], children[7], )


def edgeProc(q0, q1):
    if q0.is_leaf() or q1.is_leaf():
        return
    
    # Recursive cals to edgeProc for edge-adjacent children spanning the common edge
    q0_children = q0.get_or_copy()
    q1_children = q1.get_or_copy()

    # if q1 is to the right of q0
    if (q0.base_coords[0] != q1.base_coords[1]):
        edge_pairs = [(q0_children[1], q1_children[0]),
                    (q0_children[3], q1_children[2]),
                    (q0_children[5], q1_children[4]),
                    (q0_children[7], q1_children[6])]
        for child1, child2 in edge_pairs:
            edgeProc(child1, child2)

        vertProc(q0_children[1], q1_children[0], q0_children[3], q1_children[2], q0_children[5], q1_children[4], q0_children[7], q1_children[6])
    else:
        edge_pairs = [(q0_children[4], q1_children[0]),
                    (q0_children[5], q1_children[1]),
                    (q0_children[6], q1_children[2]),
                    (q0_children[7], q1_children[3])]
        for child1, child2 in edge_pairs:
            edgeProc(child1, child2)

        vertProc(q0_children[4], q0_children[5], q0_children[6], q0_children[7], q1_children[0], q1_children[1], q1_children[2], q1_children[3])


def vertProc(q0, q1, q2 ,q3 ,q4 ,q5 ,q6 ,q7):
    if any(q.is_leaf() for q in [q0, q1 ,q2 ,q3 ,q4 ,q5 ,q6 ,q7]):
        return
    
    # Recursive calls to vertProc for the 8 children meeting at the common
    #  vertex
    q0_children = q0.get_or_copy()
    q1_children = q1.get_or_copy()
    q2_children = q2.get_or_copy()
    q3_children = q3.get_or_copy()
    
    q4_children = q4.get_or_copy()
    q5_children = q5.get_or_copy()
    q6_children = q6.get_or_copy()
    q7_children = q7.get_or_copy()

    vertProc(q0_children[7], q1_children[6], q2_children[5], q3_children[4], q4_children[3], q5_children[2], q6_children[1], q7_children[0])


class DMCOutput:
    def __init__(self) -> None:
        self.vertices = None
        self.triangles =None
        self.vertices_len = 0

    def insert(self, verts, fcs):
        if (self.vertices is None):
            self.vertices = np.array(verts)
            self.triangles =np.array(fcs) 
            self.vertices_len = verts.shape[0]
        else:
            self.vertices = np.concatenate((self.vertices, verts), axis=0)
            fcs += self.vertices_len
            self.vertices_len += verts.shape[0]
            self.triangles = np.concatenate((self.triangles, fcs), axis=0)


def traverse(node):
    if node is None:
        return
    
    # check if the current node is a leaf
    if all(child is None for child in node.children):
        c000 = compute_f(field_values, node.base_coords[0], node.base_coords[0], node.base_coords[0])
        c001 = compute_f(field_values, node.base_coords[0], node.base_coords[0], node.base_coords[1])
        c010 = compute_f(field_values, node.base_coords[0], node.base_coords[1], node.base_coords[0])
        c011 = compute_f(field_values, node.base_coords[0], node.base_coords[1], node.base_coords[1])
        c100 = compute_f(field_values, node.base_coords[1], node.base_coords[0], node.base_coords[0])
        c101 = compute_f(field_values, node.base_coords[1], node.base_coords[0], node.base_coords[1])
        c110 = compute_f(field_values, node.base_coords[1], node.base_coords[1], node.base_coords[0])
        c111 = compute_f(field_values, node.base_coords[1], node.base_coords[1], node.base_coords[1])
        
        scalar_field = np.array([[[c000, c001],
                         [c010, c011]],
                         [[c100, c101],
                         [c110, c111]]])
        vertices, triangles = mcubes.marching_cubes(scalar_field, node.isolevel)
        # print(f'vertices:{vertices}, triangles:{triangles}')
        if(triangles.size != 0):
            output.insert(vertices, triangles)

    for child in node.children:
        traverse(child)


def dualMarchingCubes(scalar_field_values, isovalue):
    global field_values
    field_values = scalar_field_values
    global isolevel
    isolevel = isovalue
    global output
    base_coords = np.array([0,0,0])
    side_length = 1
    root = OctreeNode(base_coords, side_length)
    faceProc(root)
    output = DMCOutput()
    traverse(root)
    return output.vertices, output.triangles

def main():
    dim=4
    myfield_values = np.array([[[random.randint(-dim-2,dim+2) for i in range(dim)] for j in range(dim)] for k in range(dim)])
    isoval = 2
    vertices, triangles = dualMarchingCubes(myfield_values, isoval)
    # print(f'vertices:{vertices}, triangles:{triangles}')
main()


# def main():
#     dim = 4
#     base_coords=np.array([0,0,0])
#     side_length=1
#     global field_values
#     field_values = np.array([[[random.randint(-dim-2,dim+2) for i in range(dim)] for j in range(dim)] for k in range(dim)]) 
#     print(f'field values:{field_values}')
#     root = OctreeNode(base_coords, side_length)

#     print(f'{root.children}')
#     faceProc(root)
#     global output
#     output = DMCOutput()
#     traverse(root)
#     print(f'vertices:{output.vertices}, triangles:{output.triangles}, vertices len:{output.vertices_len}')
#     # print(f'vertices shape:{output.vertices.shape}, faces shape: {output.faces.shape}')
#     return output.vertices, output.triangles

# main()
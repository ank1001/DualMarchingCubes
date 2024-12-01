import numpy as np
# from OctreeNode4 import field_values
from LookupTable2 import tri_table
from FeatureIsolation4 import compute_f
# Edge table to map edges to vertex pairs
edge_table = [
    [0, 1], [1, 2], [2, 3], [3, 0],
    [4, 5], [5, 6], [6, 7], [7, 4],
    [0, 4], [1, 5], [2, 6], [3, 7]
]

def marching_cubes(cube_positions, cube_isovalues, isolevel):
    # Function to get the cube index
    def get_cube_index(isovalues, isolevel):
        index = 0
        for i, value in enumerate(isovalues):
            if value < isolevel:
                index |= 1 << i
        return index

    # Function to interpolate vertices
    def interpolate_vertex(p1, p2, val1, val2, isolevel):
        if abs(isolevel - val1) < 0.00001:
            return p1
        if abs(isolevel - val2) < 0.00001:
            return p2
        if abs(val1 - val2) < 0.00001:
            return p1
        t = (isolevel - val1) / (val2 - val1)
        return p1 + t * (p2 - p1)

    vertices = []
    triangles = []
    vert_dict = {}
    vert_count = 0

    # Get the cube index
    cube_index = get_cube_index(cube_isovalues, isolevel)
    print(f'cubeindex:{cube_index}')

    # Get the corresponding triangles from the table
    edges = tri_table[cube_index]
    print(f'edges:{edges}')
    if not edges:
        return np.array(vertices), np.array(triangles)
    
    # Interpolate vertices
    for i in range(0, len(edges), 3):
        e1 = edge_table[edges[i]]
        e2 = edge_table[edges[i+1]]
        e3 = edge_table[edges[i+2]]
        
        v1 = interpolate_vertex(cube_positions[e1[0]], cube_positions[e1[1]], cube_isovalues[e1[0]], cube_isovalues[e1[1]], isolevel)
        v2 = interpolate_vertex(cube_positions[e2[0]], cube_positions[e2[1]], cube_isovalues[e2[0]], cube_isovalues[e2[1]], isolevel)
        v3 = interpolate_vertex(cube_positions[e3[0]], cube_positions[e3[1]], cube_isovalues[e3[0]], cube_isovalues[e3[1]], isolevel)
        
        for v in [v1, v2, v3]:
            v_tuple = tuple(v)
            if v_tuple not in vert_dict:
                vert_dict[v_tuple] = vert_count
                vertices.append(v.tolist())
                vert_count += 1
        
        triangles.append([vert_dict[tuple(v1)], vert_dict[tuple(v2)], vert_dict[tuple(v3)]])
    
    return np.array(vertices), np.array(triangles)
# Example usage
def mod_marching_cubes(cell_params, isolevel, field_values):
    vertices_coords = np.array(cell_params[:, 1:])
    scalar_field = np.array([compute_f(field_values, p[0], p[1], p[2]) for p in vertices_coords])
    vertices, triangles = marching_cubes(vertices_coords,scalar_field, isolevel)
    return np.array(vertices), np.array(triangles)


# scalar_field = np.random.rand(3, 3, 3)
# isolevel = 0.5

# vertices, triangles = marching_cubes(scalar_field, isolevel)

# print("Vertices:", vertices)
# print("Triangles:", triangles)

# def main():
#     cell_params = np.array([[0.69869749, 0.37584364, 0.40615453, 0.28601576],
#        [0.60851223, 0.03163737, 0.1016379 , 0.2606629 ],
#        [0.70948277, 0.17368685, 0.5189909 , 0.04330313],
#        [0.7617489 , 0.11840514, 0.68552033, 0.8281011 ],
#        [0.14304385, 0.79745545, 0.6676333 , 0.91678687],
#        [0.60626532, 0.72849506, 0.10618249, 0.32032418],
#        [0.79292376, 0.86369677, 0.85736879, 0.89069221],
#        [0.83297741, 0.4932701 , 0.78915049, 0.17245434]])
    
#     isolevel = 0.15
#     vertices, faces = mod_marching_cubes(cell_params, isolevel)
#     print(vertices, faces)
# if __name__=="__main__":
#     main()
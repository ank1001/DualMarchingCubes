import numpy as np
import math
from scipy.optimize import minimize


class Feature_Isolation:

    def __init__(self, base_coords, side_length, field_values, isolevel, number_of_points) -> None:
        self.base_coords = base_coords
        self.side_length = side_length
        self.number_of_points = number_of_points
        self.res = self.compute_QEF(field_values, isolevel)
        self.QEF = self.res[0]
        self.cell_params = self.res[1]

    def compute_QEF(self, field_values, isolevel):
        mid_point = np.array([
        self.base_coords[0] + self.side_length / 2,
        self.base_coords[1] + self.side_length / 2,
        self.base_coords[2] + self.side_length / 2
        ])

        # Initial guess for optimization
        initial_guess = np.hstack([isolevel, mid_point])  # w = 0, (x, y, z) = midpoint

        # Bounds for optimization
        bounds = [
        (None, None),  # No bounds for w
        (self.base_coords[0], self.base_coords[0] + self.side_length),  # x bounds
        (self.base_coords[1], self.base_coords[1] + self.side_length),  # y bounds
        (self.base_coords[2], self.base_coords[2] + self.side_length)   # z bounds
        ]

        # Minimize the QEF
        result = minimize(
        self.compute_E,
        initial_guess,
        args=(field_values, isolevel, mid_point),
        bounds=bounds
        )

        # Return the minimized QEF and the optimal parameters
        return result.fun, result.x

    def compute_E(self, params, field_values, isolevel, mid_point):
        

        points = np.random.rand(self.number_of_points, 3)

        w, x, y, z = params

        # Compute the gradient at each point
        gradient = np.array([compute_gradient(field_values, p[0], p[1], p[2]) for p in points])

        # Compute the magnitude of each gradient vector
        gr_magnitude = np.array([np.sqrt(p[0]**2 +  p[1]**2+ p[2]**2) for p in gradient])

        # Compute the error term E for each point
        E_arr = np.array([
        ((w - np.dot(gradient[i], np.array([x - points[i][0], y - points[i][1], z - points[i][2]]))) ** 2) / (1 + gr_magnitude[i]**2)
        for i in range(len(gradient))
        ])

        # Sum up all error terms
        E = np.sum(E_arr)
        return E
    

# computes gradient at any point
def compute_gradient(field_values, x, y, z, h = 1e-5):
    dfdx = (compute_f(field_values, x + h, y, z) - compute_f(field_values, x - h, y, z)) / (2 * h)
    dfdy = (compute_f(field_values, x, y + h, z) - compute_f(field_values, x, y - h, z)) / (2 * h)
    dfdz = (compute_f(field_values, x, y, z + h) - compute_f(field_values, x, y, z - h)) / (2 * h)
    return np.array([dfdx, dfdy, dfdz])

# using trilinear interpolation
def compute_f(field_values, sx, sy, sz):
    # (x,y,z) is a point on a unit cube that will be mapped to field_values

    # if out out bounds
    if ((sx>1) | (sy>1) | (sz>1) | (sx<0) | (sy<0) | (sz<0)):
        interpolated_value= 0
    else:
        dimx, dimy, dimz = field_values.shape

        # base coords of cube
        x, y, z = sx*(dimx-1), sy*(dimy-1), sz*(dimz-1)
        

        i0, j0, k0 = math.floor(x), math.floor(y), math.floor(z)
        if(i0==dimx-1):
            i0 = i0-1
        if(j0==dimy-1):
            j0 = j0-1
        if(k0==dimz-1):
            k0 = k0-1
        i1, j1, k1 = i0+1, j0+1, k0+1
        wx = (x - i0) / (i1 - i0) if i1 != i0 else 0
        wy = (y - j0) / (j1 - j0) if j1 != j0 else 0
        wz = (z - k0) / (k1 - k0) if k1 != k0 else 0

        f000 = field_values[i0, j0, k0]
        f100 = field_values[i1, j0, k0]
        f010 = field_values[i0, j1, k0]
        f110 = field_values[i1, j1, k0]
        f001 = field_values[i0, j0, k1]
        f101 = field_values[i1, j0, k1]
        f011 = field_values[i0, j1, k1]
        f111 = field_values[i1, j1, k1]

        interpolated_value = (
        (1 - wx) * (1 - wy) * (1 - wz) * f000 +
        wx * (1 - wy) * (1 - wz) * f100 +
        (1 - wx) * wy * (1 - wz) * f010 +
        wx * wy * (1 - wz) * f110 +
        (1 - wx) * (1 - wy) * wz * f001 +
        wx * (1 - wy) * wz * f101 +
        (1 - wx) * wy * wz * f011 +
        wx * wy * wz * f111
        )
    return interpolated_value


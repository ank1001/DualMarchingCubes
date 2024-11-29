import numpy as np
import math
from scipy.optimize import minimize


class Feature_Isolation:

    def __init__(self, base_coords, side_length, field_values, isolevel) -> None:
        self.base_coords = base_coords
        self.side_length = side_length
        self.res = self.compute_QEF(field_values, isolevel)
        self.QEF = self.res[0]
        self.grid_index = self.res[1]
        self.QEF_threshold = 0 # each cell must have QEF value less than -100

    def compute_QEF(self, field_values, isolevel):
        mid_point = np.array([self.base_coords[0] + self.side_length/2, self.base_coords[1] + self.side_length/2,self.base_coords[2] + self.side_length/2])
        initial_guess = mid_point
        bounds = [(self.base_coords[0], self.base_coords[0]+self.side_length), (self.base_coords[1], self.base_coords[1]+self.side_length), (self.base_coords[2], self.base_coords[2]+self.side_length)]
        # print(f'bounds:{bounds}')
        result = minimize(self.compute_E, initial_guess, args=(field_values, isolevel, mid_point), bounds = bounds)
        return (result.fun, result.x)

    def compute_E(self, params, field_values, isolevel, mid_point):
        # mid point 
        xi, yi, zi = mid_point
        x,y,z = params
        gradient = compute_gradient(field_values, xi, yi, zi)
        gr_magnitude = np.sqrt(gradient[0]**2 + gradient[1]**2 + gradient[2]**2)
        E = ((isolevel - ( np.dot(gradient, np.array([x-xi, y-yi, z-zi])) ))**2)/(1 + gr_magnitude**2 )
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



    # # using trilinear interpolation
    # def compute_f2(self, field_values, sx, sy, sz):
    #     # (x,y,z) is a point on a unit cube that will be mapped to field_values

    #     # if out out bounds
    #     if ((sx>1) | (sy>1) | (sz>1) | (sx<0) | (sy<0) | (sz<0)):
    #         interpolated_value= 0
    #     else:
    #         dimx, dimy, dimz = field_values.shape

    #         # base coords of cube
    #         x, y, z = sx*(dimx-1), sy*(dimy-1), sz*(dimz-1)
            

    #         i0, j0, k0 = math.floor(x), math.floor(y), math.floor(z)
    #         if(i0==dimx-1):
    #             i0 = i0-1
    #         if(j0==dimy-1):
    #             j0 = j0-1
    #         if(k0==dimz-1):
    #             k0 = k0-1
    #         i1, j1, k1 = i0+1, j0+1, k0+1
    #         id, jd, kd = (x-i0)/(i1-i0), (y-j0)/(j1-j0), (z-k0)/(k1-k0)

    #         c000 = field_values[i0, j0, k0]
    #         c001 = field_values[i0, j0, k1]
    #         c010 = field_values[i0, j1, k0]
    #         c011 = field_values[i0, j1, k1]
    #         c100 = field_values[i1, j0, k0]
    #         c101 = field_values[i1, j0, k1]
    #         c110 = field_values[i1, j1, k0]
    #         c111 = field_values[i1, j1, k1]

    #         c00 = c000*(1-id) + c100*id
    #         c01 = c001*(1-id) + c101*id
    #         c10 = c010*(1-id) + c110*id
    #         c11 = c011*(1-id) + c111*id

    #         c0 = c00*(1-jd) + c10*jd
    #         c1 = c01*(1-jd) + c11*jd

    #         interpolated_value = c0*(1-kd) + c1*kd
    #     return interpolated_value


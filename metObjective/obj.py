
import numpy as np

# Flux parameters
TARGET_FLUX1 = 1 # target v for first flux
TARGET_FLUX2 = 10 # target f for second flux

TARGET_TOTAL_FLUX = 10 # constraint on total flux (sum of first and second)

# Optimization parameters
MAX_ITERATIONS = 1000
INITIAL_STEPSIZE = 1e-3
STEP_SUCCESS_FOLDCHANGE = 1.1
STEP_FAILURE_FOLDCHANGE = 0.5
GRADIENT_PRECISION = 1e-3

def nullspace(matrix):
    # Finds the nullspace of a matrix using the SVD
    (s, u, vt) = np.linalg.svd(matrix, full_matrices = True)
    iszero = np.where(u == 0)

    if np.any(iszero):
        startat = np.where(iszero)

    else:
        startat = u.size

    return vt[:, startat:]

def add_to_coordinate(vector, i, dx):
    # Adds a value to a specific position in a vector
    out = vector.copy()

    out[i] += dx

    return out

def approx_gradient(f, x, r):
    # Compute an approximate gradient using a symmetric difference
    return np.array([
        (f(add_to_coordinate(x, i, +r)) - f(add_to_coordinate(x, i, -r)))/(2*r)
        for i in range(x.size)
        ])

# Define objective functions

def abs_residual(current, expected):
    return current - expected

def rel_residual(current, expected):
    return 1 - current/expected

def log_residual(current, expected):
    return np.log(current / expected)

def l1_norm(residuals):
    return np.sum(np.abs(residuals))

def l2_norm(residuals): # actually the squared L2 norm
    return np.sum(np.square(residuals))

# Convert basis variables (from 'x' to 'y') to eliminate equality constraint

A = np.array([[1, 1]], np.float64)
b = np.array([TARGET_TOTAL_FLUX], np.float64)

x0 = np.linalg.lstsq(A, b, rcond = -1)[0]

N = nullspace(A)

y0 = np.zeros(N.shape[1], np.float64)

expected_flux = np.array([TARGET_FLUX1, TARGET_FLUX2], np.float64)

# Convenience functions for switching between x and y
y_to_x = lambda y: N.dot(y) + x0
x_to_y = lambda x: np.linalg.lstsq(N, x - x0, rcond = None)[0]

# Solve each optimization problem

for residual in (abs_residual, rel_residual, log_residual):
    for norm in (l1_norm, l2_norm):
        name = '{}, {}'.format(residual.__name__, norm.__name__)
        objective_function = lambda x: norm(residual(x, expected_flux))

        y = y0.copy()
        stepsize = INITIAL_STEPSIZE
        f = objective_function(y_to_x(y))

        for iterate in range(MAX_ITERATIONS):
            grad = x_to_y(approx_gradient(objective_function, y_to_x(y), GRADIENT_PRECISION))

            y_new = y - stepsize * grad
            f_new = objective_function(y_to_x(y_new))

            if np.isfinite(f_new) and (f_new < f):
                y = y_new
                f = f_new

                stepsize *= STEP_SUCCESS_FOLDCHANGE

            else:
                stepsize *= STEP_FAILURE_FOLDCHANGE

        print('{}: v1 = {:0.2f}, v2 = {:0.2f}'.format(name, *y_to_x(y)))

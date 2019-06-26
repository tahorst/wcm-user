'''
Explores solution of transcription unit problem for different cases of gene to
TU mapping.

1. perfectly determined: n TU = m genes, full rank A matrix - direct solution
2. under determined: n TU > m genes, fat A matrix - least norm solution
3. over determined: n TU < m genes, skinny A matrix - least squares solution
'''

import numpy as np

# Square A matrix
# Full rank, determined
A1 = np.array([
	[1, 0],
	[1, 1]])
# Full rank, under determined
A2 = np.array([
	[1, 1, 0],
	[1, 0, 1]])
# Full rank, over determined
A3 = np.array([
	[1],
	[1]])

A = np.zeros((6, 6))
A[:2, :2] = A1
A[2:4, 2:5] = A2
A[4:, 5:] = A3

e1 = np.array([1, 2])
e2 = np.array([3, 2])
# e2 = np.array([3, 0.5])  # to get negative from least norm
e3 = np.array([2.5, 3.5])
exp = np.hstack((e1, e2, e3))

# Direct solution, requires sub matrices to be full rank
tu1 = np.linalg.inv(A1).dot(e1)
tu2 = A2.T.dot(np.linalg.inv(A2.dot(A2.T))).dot(e2)
tu3 = np.linalg.inv(A3.T.dot(A3)).dot(A3.T).dot(e3)
tu = np.hstack((tu1, tu2, tu3))

# Other solutions, A is not full rank - give same results as tu
tu_pinv = np.linalg.pinv(A).dot(exp)
tu_lstsq, _, _, _ = np.linalg.lstsq(A, exp, rcond=None)

print('\nSquare A:\n{}'.format(A))
print('submatrix solution:\n{}'.format(tu))
print('pinv:\n{}'.format(tu_pinv))
print('lstsq:\n{}'.format(tu_lstsq))

print('\nOriginal exp:\n{}'.format(exp))
print('Unit exp from solution:\n{}'.format(A.dot(tu)))

# Fat matrix for combined A
# Full rank, under determined
A2 = np.array([
	[1, 1, 0, 1, 0],
	[1, 0, 1, 0, 0],
	[1, 1, 1, 0, 1]])
# e2 = np.array([3, 2, 2.5])  # to get negative from least norm
e2 = np.array([3, 2, 3.5])
A = np.zeros((7, 8))
A[:2, :2] = A1
A[2:5, 2:7] = A2
A[5:, 7:] = A3
exp = np.hstack((e1, e2, e3))

# Direct solution, requires sub matrices to be full rank
tu1 = np.linalg.inv(A1).dot(e1)
tu2 = A2.T.dot(np.linalg.inv(A2.dot(A2.T))).dot(e2)
tu3 = np.linalg.inv(A3.T.dot(A3)).dot(A3.T).dot(e3)
tu = np.hstack((tu1, tu2, tu3))

# Other solutions, A is not full rank - give same results as tu
tu_pinv = np.linalg.pinv(A).dot(exp)
tu_lstsq, _, _, _ = np.linalg.lstsq(A, exp, rcond=None)

print('\nFat A:\n{}'.format(A))
print('submatrix solution:\n{}'.format(tu))
print('pinv:\n{}'.format(tu_pinv))
print('lstsq:\n{}'.format(tu_lstsq))

print('\nOriginal exp:\n{}'.format(exp))
print('Unit exp from solution:\n{}'.format(A.dot(tu)))

import ipdb; ipdb.set_trace()

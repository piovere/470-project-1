from numpy.linalg import inv, solve
from numpy import zeros, dot, array, reshape


sigma = 0.15
W = 10.0
S = 10 ** 8
num_nodes = 10000

def set_matrix(nodes, geometry, D, w=W, sigma=sigma):
	deltaX = w / (nodes - 1.0)
	n = nodes-1

	matrix = zeros((n, n))
	if geometry == "cartesian":
		matrix[0][0] = D / deltaX ** 2.0 + 0.5 * sigma
		matrix[0][1] = -1.0 * D / deltaX ** 2.0
		for i in range(1, n-1):
			matrix[i][i] = 2 * D / deltaX ** 2 + sigma
			matrix[i][i-1] = -1.0 * D / deltaX ** 2
			matrix[i][i+1] = -1.0 * D / deltaX ** 2
		matrix[n-1][n-1] = 2 * D / deltaX ** 2.0 + sigma
		matrix[n-1][n-2] = -1.0 * D / deltaX ** 2.0

	return matrix

A = set_matrix(num_nodes, "cartesian", 9.0)

# print("This is the matrix")
# print(A)

# print("This is the inverted matrix")
# print(inv(A))

# print("Verify that they are inverses")
# print(dot(A, inv(A)))

n = num_nodes - 1
S_vector = zeros(n)
S_vector[0] = S / (2.0 * W / n)

# print("This is the S vector")
# print(S_vector)

print("Multiplying S vector by inverse A")
solution = dot(inv(A), S_vector)
print(solution[0])
print(solution[n-1])

# print("What happens if numpy solves it for me?")
# phi = solve(A, S_vector)
# print(phi)

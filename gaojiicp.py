import numpy as np
from scipy.spatial import KDTree

DISTANCE_THRESHOLD = 1.0  # 设定一个距离阈值来剔除噪声点

def compute_centroid(points):
    return np.mean(points, axis=0)

def compute_icp(A, B, max_iterations=100, tolerance=1e-6):
    tree = KDTree(B)
    R = np.eye(3)
    t = np.zeros(3)
    
    for iteration in range(max_iterations):
        # Apply current transformation to A
        A_transformed = (R @ A.T).T + t
        
        # Find nearest neighbors in B for each point in A
        distances, indices = tree.query(A_transformed)
        
        # Filter out pairs with distance greater than the threshold
        valid_pairs = distances < DISTANCE_THRESHOLD
        A_valid = A[valid_pairs]
        B_valid = B[indices[valid_pairs]]
        
        if len(A_valid) < 3:
            print("Not enough valid pairs")
            break
        
        # Compute centroids of valid pairs
        centroid_A = compute_centroid(A_valid)
        centroid_B = compute_centroid(B_valid)
        
        # Center the points
        A_centered = A_valid - centroid_A
        B_centered = B_valid - centroid_B
        
        # Compute covariance matrix
        H = A_centered.T @ B_centered
        
        # Singular Value Decomposition
        U, S, Vt = np.linalg.svd(H)
        
        # Compute rotation
        R_new = Vt.T @ U.T
        if np.linalg.det(R_new) < 0:
            Vt[2, :] *= -1
            R_new = Vt.T @ U.T
        
        # Compute translation
        t_new = centroid_B - R_new @ centroid_A
        
        # Check for convergence
        if np.linalg.norm(R_new - R) < tolerance and np.linalg.norm(t_new - t) < tolerance:
            break
        
        R = R_new
        t = t_new
    
    return R, t

# 示例点云数据
A = np.array([
    [1, 2, 3], [4, 5, 6], [7, 8, 9],
    [2, 3, 4], [5, 6, 7], [8, 9, 10],
    [3, 4, 5], [6, 7, 8], [9, 10, 11], [10, 11, 12]
])

B = np.array([
    [1.2, 2.1, 2.9], [4.1, 5.0, 6.2], [6.9, 8.2, 9.1],
    [2.2, 3.1, 4.3], [5.1, 6.0, 7.2], [8.2, 9.1, 10.0],
    [3.1, 4.0, 5.3], [6.0, 7.2, 8.1], [9.2, 10.1, 10.9],
    [10.0, 11.1, 12.0], [13.0, 14.1, 15.2], [16.0, 17.1, 18.0]
])

# 执行ICP算法
R, t = compute_icp(A, B)

print("旋转矩阵 R:")
print(R)
print("平移向量 t:")
print(t)
#------------------------
import numpy as np

class KDTree:
    def __init__(self, points):
        self.k = points.shape[1]
        self.root = self.build_tree(points)

    class Node:
        def __init__(self, point, left, right):
            self.point = point
            self.left = left
            self.right = right

    def build_tree(self, points, depth=0):
        if len(points) == 0:
            return None
        
        axis = depth % self.k
        points = points[points[:, axis].argsort()]
        median = len(points) // 2

        return self.Node(
            point=points[median],
            left=self.build_tree(points[:median], depth + 1),
            right=self.build_tree(points[median + 1:], depth + 1)
        )

    def nearest(self, point, return_distance=True):
        best = [None, float('inf')]
        def recursive_search(node, depth=0):
            if node is None:
                return

            axis = depth % self.k
            dist = np.linalg.norm(point - node.point)

            if dist < best[1]:
                best[0] = node.point
                best[1] = dist

            diff = point[axis] - node.point[axis]
            close, away = (node.left, node.right) if diff < 0 else (node.right, node.left)

            recursive_search(close, depth + 1)
            if diff ** 2 < best[1]:
                recursive_search(away, depth + 1)

        recursive_search(self.root)
        return best if return_distance else best[0]

#### ICP算法实现

```python
import numpy as np

DISTANCE_THRESHOLD = 1.0

class KDTree:
    def __init__(self, points):
        self.k = points.shape[1]
        self.root = self.build_tree(points)

    class Node:
        def __init__(self, point, left, right):
            self.point = point
            self.left = left
            self.right = right

    def build_tree(self, points, depth=0):
        if len(points) == 0:
            return None
        
        axis = depth % self.k
        points = points[points[:, axis].argsort()]
        median = len(points) // 2

        return self.Node(
            point=points[median],
            left=self.build_tree(points[:median], depth + 1),
            right=self.build_tree(points[median + 1:], depth + 1)
        )

    def nearest(self, point, return_distance=True):
        best = [None, float('inf')]
        def recursive_search(node, depth=0):
            if node is None:
                return

            axis = depth % self.k
            dist = np.linalg.norm(point - node.point)

            if dist < best[1]:
                best[0] = node.point
                best[1] = dist

            diff = point[axis] - node.point[axis]
            close, away = (node.left, node.right) if diff < 0 else (node.right, node.left)

            recursive_search(close, depth + 1)
            if diff ** 2 < best[1]:
                recursive_search(away, depth + 1)

        recursive_search(self.root)
        return best if return_distance else best[0]

def compute_centroid(points):
    return np.mean(points, axis=0)

def compute_icp(A, B, max_iterations=100, tolerance=1e-6):
    tree = KDTree(B)
    R = np.eye(3)
    t = np.zeros(3)
    
    for iteration in range(max_iterations):
        A_transformed = (R @ A.T).T + t
        valid_pairs = []
        
        for point in A_transformed:
            nearest_point, distance = tree.nearest(point)
            if distance < DISTANCE_THRESHOLD:
                valid_pairs.append((point, nearest_point))
        
        if len(valid_pairs) < 3:
            print("Not enough valid pairs")
            break
        
        A_valid = np.array([p[0] for p in valid_pairs])
        B_valid = np.array([p[1] for p in valid_pairs])
        
        centroid_A = compute_centroid(A_valid)
        centroid_B = compute_centroid(B_valid)
        
        A_centered = A_valid - centroid_A
        B_centered = B_valid - centroid_B
        
        H = A_centered.T @ B_centered
        U, S, Vt = np.linalg.svd(H)
        
        R_new = Vt.T @ U.T
        if np.linalg.det(R_new) < 0:
            Vt[2, :] *= -1
            R_new = Vt.T @ U.T
        
        t_new = centroid_B - R_new @ centroid_A
        
        if np.linalg.norm(R_new - R) < tolerance and np.linalg.norm(t_new - t) < tolerance:
            break
        
        R = R_new
        t = t_new
    
    return R, t

# 示例数据
A = np.array([
    [1, 2, 3], [4, 5, 6], [7, 8, 9],
    [2, 3, 4], [5, 6, 7], [8, 9, 10],
    [3, 4, 5], [6, 7, 8], [9, 10, 11], [10, 11, 12]
])

B = np.array([
    [1.2, 2.1, 2.9], [4.1, 5.0, 6.2], [6.9, 8.2, 9.1],
    [2.2, 3.1, 4.3], [5.1, 6.0, 7.2], [8.2, 9.1, 10.0],
    [3.1, 4.0, 5.3], [6.0, 7.2, 8.1], [9.2, 10.1, 10.9],
    [10.0, 11.1, 12.0], [13.0, 14.1, 15.2], [16.0, 17.1, 18.0]
])

# 执行ICP算法
R, t = compute_icp(A, B)

print("旋转矩阵 R:")
print(R)
print("平移向量 t:")
print(t)

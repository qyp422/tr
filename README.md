KD-Tree实现:

使用静态数组和递归方法构建KD-Tree。
nearest函数用于找到最近的点。
ICP算法:

构建KD-Tree。
使用迭代最近点算法找到每个点的最近邻点，并计算质心和旋转矩阵。
辅助函数:

computeCentroid计算质心。
subtractCentroid从点中减去质心。
computeSVD计算奇异值分解。
matrixMultiply和matrixTranspose进行矩阵运算。
matrixVectorMultiply进行矩阵与向量的乘法。
这个代码实现了一个基本的ICP算法，能够处理带有噪声的点云数据。通过自定义的KD-Tree实现，提供了高效的最近邻搜索功能。







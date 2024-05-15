#include <iostream>
#include <cmath>

const int MAX_POINTS = 1000;
const int K = 3; // KD 树中每个节点的维度

// 点的类定义
class Point {
public:
    int id;
    double coords[K]; // 存储点的坐标
    int cluster; // 存储点的聚类标签

    Point(int _id, double x, double y, double z) : id(_id), cluster(-1) {
        coords[0] = x;
        coords[1] = y;
        coords[2] = z;
    }
};

// KD 树节点定义
class KDNode {
public:
    Point* point;
    KDNode* left;
    KDNode* right;

    KDNode(Point* p) : point(p), left(nullptr), right(nullptr) {}
};

// KD 树类定义
class KDTree {
public:
    KDTree() : root(nullptr) {}

    // 添加一个点到 KD 树中
    void insert(Point* point) {
        root = insertRec(root, point, 0);
    }

    // KD 树搜索
    void search(double cutoff, Point* points[MAX_POINTS]) {
        searchRec(root, cutoff, points);
    }

private:
    KDNode* root;

    // 在 KD 树中递归插入一个点
    KDNode* insertRec(KDNode* node, Point* point, int depth) {
        if (node == nullptr) {
            return new KDNode(point);
        }

        int axis = depth % K;
        if (point->coords[axis] < node->point->coords[axis]) {
            node->left = insertRec(node->left, point, depth + 1);
        } else {
            node->right = insertRec(node->right, point, depth + 1);
        }

        return node;
    }

    // 在 KD 树中搜索并进行聚类
    void searchRec(KDNode* node, double cutoff, Point* points[MAX_POINTS]) {
        if (node == nullptr) {
            return;
        }

        // 遍历节点的所有邻居，如果距离小于截断距离，则将其分配到同一类别
        for (int i = 0; i < MAX_POINTS; ++i) {
            if (points[i] != nullptr && points[i]->cluster == -1 && distance(node->point, points[i]) < cutoff) {
                points[i]->cluster = points[i]->id;
            }
        }

        // 递归搜索左子树和右子树
        searchRec(node->left, cutoff, points);
        searchRec(node->right, cutoff, points);
    }

    // 计算两个点之间的空间距离
    double distance(Point* p1, Point* p2) {
        double dist = 0.0;
        for (int i = 0; i < K; ++i) {
            dist += std::pow(p1->coords[i] - p2->coords[i], 2);
        }
        return std::sqrt(dist);
    }
};

int main() {
    // 创建 KD 树
    KDTree tree;

    // 添加一些点到 KD 树中
    tree.insert(new Point(1, 0.0, 0.0, 0.0));
    tree.insert(new Point(2, 1.0, 1.0, 1.0));
    tree.insert(new Point(3, 2.0, 2.0, 2.0));
    tree.insert(new Point(4, 3.0, 3.0, 3.0));
    tree.insert(new Point(5, 5.0, 5.0, 5.0));

    // 创建点数组用于存储聚类结果
    Point* points[MAX_POINTS] = {nullptr};
    points[0] = new Point(1, 0.0, 0.0, 0.0);
    points[1] = new Point(2, 1.0, 1.0, 1.0);
    points[2] = new Point(3, 2.0, 2.0, 2.0);
    points[3] = new Point(4, 3.0, 3.0, 3.0);
    points[4] = new Point(5, 5.0, 5.0, 5.0);

    // 执行 KD 树搜索并进行聚类
    tree.search(2.0, points);

    // 输出聚类结果
    for (int i = 0; i < MAX_POINTS; ++i) {
        if (points[i] != nullptr) {
            std::cout << "Point " << points[i]->id << " belongs to cluster " << points[i]->cluster << std::endl;
        }
    }

    return 0;
}

//-----------------------------------------------------------
#include <iostream>
#include <cmath>

const int MAX_POINTS = 1000;
const int K = 3; // KD 树中每个节点的维度

// 点的类定义
class Point {
public:
    int id;
    double coords[K]; // 存储点的坐标

    Point(int _id, double x, double y, double z) : id(_id) {
        coords[0] = x;
        coords[1] = y;
        coords[2] = z;
    }
};

// KD 树节点定义
class KDNode {
public:
    Point* point;
    KDNode* left;
    KDNode* right;

    KDNode(Point* p) : point(p), left(nullptr), right(nullptr) {}
};

// KD 树类定义
class KDTree {
public:
    KDTree() : root(nullptr) {}

    // 添加一个点到 KD 树中
    void insert(Point* point) {
        root = insertRec(root, point, 0);
    }

    // 执行 KD 树搜索并进行聚类
    void searchAndCluster(double cutoff) {
        int currentLabel = 0;
        searchAndClusterRec(root, cutoff, currentLabel);
    }

private:
    KDNode* root;

    // 在 KD 树中递归插入一个点
    KDNode* insertRec(KDNode* node, Point* point, int depth) {
        if (node == nullptr) {
            return new KDNode(point);
        }

        int axis = depth % K;
        if (point->coords[axis] < node->point->coords[axis]) {
            node->left = insertRec(node->left, point, depth + 1);
        } else {
            node->right = insertRec(node->right, point, depth + 1);
        }

        return node;
    }

    // 在 KD 树中搜索并进行聚类
    void searchAndClusterRec(KDNode* node, double cutoff, int& currentLabel) {
        if (node == nullptr) {
            return;
        }

        // 遍历节点的所有邻居，如果距离小于截断距离，则将其分配到同一类别
        for (int i = 0; i < MAX_POINTS; ++i) {
            if (points[i] != nullptr && points[i]->cluster == -1 && distance(node->point, points[i]) < cutoff) {
                points[i]->cluster = currentLabel;
            }
        }

        // 递归搜索左子树和右子树
        searchAndClusterRec(node->left, cutoff, currentLabel);
        searchAndClusterRec(node->right, cutoff, currentLabel);
    }

    // 计算两个点之间的空间距离
    double distance(Point* p1, Point* p2) {
        double dist = 0.0;
        for (int i = 0; i < K; ++i) {
            dist += std::pow(p1->coords[i] - p2->coords[i], 2);
        }
        return std::sqrt(dist);
    }
};

int main() {
    // 创建 KD 树
    KDTree tree;

    // 添加一些点到 KD 树中
    tree.insert(new Point(1, 0.0, 0.0, 0.0));
    tree.insert(new Point(2, 1.0, 1.0, 1.0));
    tree.insert(new Point(3, 2.0, 2.0, 2.0));
    tree.insert(new Point(4, 3.0, 3.0, 3.0));
    tree.insert(new Point(5, 5.0, 5.0, 5.0));

    // 执行 KD 树搜索并进行聚类
    tree.searchAndCluster(2.0);

    return 0;
}
//-----------------------------------------------------------------

#include <iostream>
#include <cmath>

const int MAX_POINTS = 1000;

// 点的类定义
class Point {
public:
    int id;
    double x, y, z, t;

    Point(int _id, double _x, double _y, double _z, double _t) : id(_id), x(_x), y(_y), z(_z), t(_t) {}
    
    // 计算点之间的空间距离
    double distance(const Point& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        double dz = z - other.z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }
};

// 聚类器类定义
class Clusterer {
public:
    // 构造函数，接受截断距离作为参数
    Clusterer(double cutoff) : cutoffDistance(cutoff), numPoints(0), currentLabel(0) {}

    // 添加一个点到聚类器中
    void addPoint(const Point& point) {
        points[numPoints++] = point;
    }

    // 执行聚类操作
    void cluster() {
        // 初始化每个点的类别为 -1，表示未分配
        int labels[MAX_POINTS];
        for (int i = 0; i < numPoints; ++i) {
            labels[i] = -1;
        }

        // 遍历每个点
        for (int i = 0; i < numPoints; ++i) {
            // 如果该点尚未分配类别
            if (labels[i] == -1) {
                // 将该点作为新的类别的代表
                labels[i] = currentLabel++;
                // 对于该点的所有邻居，如果距离小于截断距离，则将其分配到同一类别
                for (int j = i + 1; j < numPoints; ++j) {
                    if (points[i].distance(points[j]) < cutoffDistance) {
                        labels[j] = labels[i];
                    }
                }
            }
        }

        // 输出结果
        for (int i = 0; i < numPoints; ++i) {
            std::cout << "Point " << points[i].id << " belongs to cluster " << labels[i] << std::endl;
        }
    }

private:
    Point points[MAX_POINTS]; // 存储所有点
    double cutoffDistance; // 截断距离
    int numPoints; // 点的数量
    int currentLabel; // 当前类别标签
};

int main() {
    // 创建聚类器，设置截断距离为 2.0
    Clusterer clusterer(2.0);

    // 添加一些点到聚类器中
    clusterer.addPoint(Point(1, 0.0, 0.0, 0.0, 0.0));
    clusterer.addPoint(Point(2, 1.0, 1.0, 1.0, 1.0));
    clusterer.addPoint(Point(3, 2.0, 2.0, 2.0, 2.0));
    clusterer.addPoint(Point(4, 3.0, 3.0, 3.0, 3.0));
    clusterer.addPoint(Point(5, 5.0, 5.0, 5.0, 5.0));

    // 执行聚类操作
    clusterer.cluster();

    return 0;
}

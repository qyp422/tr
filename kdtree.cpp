#include <iostream>
#include <algorithm>
#include <cmath>
#include <climits>

#define MAX_POINTS 1000

struct Point {
    double x, y, z;
};

struct Node {
    Point pt; // 节点的点
    Node* left;
    Node* right;
};
// class
class KDTree {
public:
    KDTree() : root(nullptr), size(0) {}
    ~KDTree() { clear(root); }

    void build(Point points[], int n);
    Point nearest(const Point& point) const;

private:
    Node* root;
    int size; // 点的数量

    Node* buildRec(Point points[], int start, int end, int depth);
    void clear(Node* node);
    Point nearestRec(Node* node, const Point& point, int depth, Point& best, double& bestDist) const;
    double squaredDist(const Point& a, const Point& b) const;
    static bool compX(const Point& a, const Point& b) { return a.x < b.x; }
    static bool compY(const Point& a, const Point& b) { return a.y < b.y; }
    static bool compZ(const Point& a, const Point& b) { return a.z < b.z; }
};

void KDTree::build(Point points[], int n) {
    this->size = n;
    root = buildRec(points, 0, n - 1, 0);
}

Node* KDTree::buildRec(Point points[], int start, int end, int depth) {
    if (start > end) return nullptr;

    int axis = depth % 3;
    int mid = start + (end - start) / 2;

    std::nth_element(points + start, points + mid, points + end + 1, [axis](const Point& a, const Point& b) {
        if (axis == 0) return a.x < b.x;
        else if (axis == 1) return a.y < b.y;
        else return a.z < b.z;
    });

    Node* node = new Node();
    node->pt = points[mid];
    node->left = buildRec(points, start, mid - 1, depth + 1);
    node->right = buildRec(points, mid + 1, end, depth + 1);

    return node;
}

Point KDTree::nearest(const Point& point) const {
    Point best = root->pt;
    double bestDist = squaredDist(point, best);
    nearestRec(root, point, 0, best, bestDist);
    return best;
}
// gouzao
Point KDTree::nearestRec(Node* node, const Point& point, int depth, Point& best, double& bestDist) const {
    if (node == nullptr) return best;

    int axis = depth % 3;
    double d = squaredDist(point, node->pt);
    double diff = (axis == 0) ? point.x - node->pt.x : (axis == 1) ? point.y - node->pt.y : point.z - node->pt.z;
    double diff2 = diff * diff;

    if (d < bestDist) {
        bestDist = d;
        best = node->pt;
    }

    Node* near = (diff > 0) ? node->left : node->right;
    Node* away = (diff > 0) ? node->right : node->left;

    nearestRec(near, point, depth + 1, best, bestDist);

    if (diff2 < bestDist) {
        nearestRec(away, point, depth + 1, best, bestDist);
    }

    return best;
}

void KDTree::clear(Node* node) {
    if (node != nullptr) {
        clear(node->left);
        clear(node->right);
        delete node;
    }
}

double KDTree::squaredDist(const Point& a, const Point& b) const {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

int main() {
    Point points[MAX_POINTS] = {
        {2.1, 3.1, 2.4}, {5.4, 2.3, 4.1}, {9.2, 9.8, 1.1}, {4.1, 1.2, 0.3}, {6.5, 4.0, 5.1},
        {3.3, 7.6, 4.2}, {7.7, 2.1, 9.3}, {8.9, 3.4, 2.2}
    };

    KDTree tree;
    tree.build(points, 8);

    Point query = {7.5, 7.5, 7.5};
    Point nearest = tree.nearest(query);

    std::cout << "Nearest point to (" << query.x << ", " << query.y << ", " << query.z << ") is ("
              << nearest.x << ", " << nearest.y << ", " << nearest.z << ")." << std::endl;

    return 0;
}

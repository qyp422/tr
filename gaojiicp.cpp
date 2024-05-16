#include <iostream>
#include <array>
#include <cmath>
#include <limits>
#include <algorithm>

const int K = 3;
const double DISTANCE_THRESHOLD = 1.0;
const size_t MAX_POINTS = 1000;

struct Point {
    std::array<double, K> coords;
};

struct Node {
    Point point;
    Node* left;
    Node* right;
};

class KDTree {
public:
    KDTree(Point points[], size_t size) {
        root = buildTree(points, size, 0);
    }

    ~KDTree() {
        destroyTree(root);
    }

    Point nearest(const Point& target, double& best_dist) const {
        best_dist = std::numeric_limits<double>::max();
        Point best_point;
        nearest(root, target, 0, best_point, best_dist);
        return best_point;
    }

private:
    Node* root;

    Node* buildTree(Point points[], size_t size, int depth) {
        if (size == 0) return nullptr;

        size_t axis = depth % K;
        std::sort(points, points + size, [axis](const Point& a, const Point& b) {
            return a.coords[axis] < b.coords[axis];
        });

        size_t median = size / 2;
        Node* node = new Node{ points[median], nullptr, nullptr };
        node->left = buildTree(points, median, depth + 1);
        node->right = buildTree(points + median + 1, size - median - 1, depth + 1);

        return node;
    }

    void destroyTree(Node* node) {
        if (node) {
            destroyTree(node->left);
            destroyTree(node->right);
            delete node;
        }
    }

    void nearest(Node* node, const Point& target, int depth, Point& best_point, double& best_dist) const {
        if (!node) return;

        double dist = distance(node->point, target);
        if (dist < best_dist) {
            best_dist = dist;
            best_point = node->point;
        }

        size_t axis = depth % K;
        double diff = target.coords[axis] - node->point.coords[axis];
        Node* near_node = diff < 0 ? node->left : node->right;
        Node* away_node = diff < 0 ? node->right : node->left;

        nearest(near_node, target, depth + 1, best_point, best_dist);
        if (diff * diff < best_dist) {
            nearest(away_node, target, depth + 1, best_point, best_dist);
        }
    }

    double distance(const Point& a, const Point& b) const {
        double dist = 0;
        for (size_t i = 0; i < K; ++i) {
            dist += (a.coords[i] - b.coords[i]) * (a.coords[i] - b.coords[i]);
        }
        return std::sqrt(dist);
    }
};

void computeCentroid(Point points[], size_t size, double centroid[K]) {
    for (size_t i = 0; i < K; ++i) {
        centroid[i] = 0;
    }
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < K; ++j) {
            centroid[j] += points[i].coords[j];
        }
    }
    for (size_t i = 0; i < K; ++i) {
        centroid[i] /= size;
    }
}

void subtractCentroid(Point points[], size_t size, const double centroid[K]) {
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < K; ++j) {
            points[i].coords[j] -= centroid[j];
        }
    }
}

void computeSVD(double H[K][K], double U[K][K], double S[K], double V[K][K]) {
    // A simple placeholder for SVD computation. For simplicity, we assume H is a diagonal matrix.
    for (size_t i = 0; i < K; ++i) {
        for (size_t j = 0; j < K; ++j) {
            U[i][j] = (i == j) ? 1 : 0;
            V[i][j] = (i == j) ? 1 : 0;
        }
        S[i] = H[i][i];
    }
}

void matrixMultiply(double A[K][K], double B[K][K], double result[K][K]) {
    for (size_t i = 0; i < K; ++i) {
        for (size_t j = 0; j < K; ++j) {
            result[i][j] = 0;
            for (size_t k = 0; k < K; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void matrixTranspose(double A[K][K], double result[K][K]) {
    for (size_t i = 0; i < K; ++i) {
        for (size_t j = 0; j < K; ++j) {
            result[j][i] = A[i][j];
        }
    }
}

void matrixVectorMultiply(double A[K][K], const double vec[K], double result[K]) {
    for (size_t i = 0; i < K; ++i) {
        result[i] = 0;
        for (size_t j = 0; j < K; ++j) {
            result[i] += A[i][j] * vec[j];
        }
    }
}

void computeICP(Point A[], Point B[], size_t sizeA, size_t sizeB, double R[K][K], double t[K], int max_iterations = 100, double tolerance = 1e-6) {
    KDTree tree(B, sizeB);
    for (size_t i = 0; i < K; ++i) {
        for (size_t j = 0; j < K; ++j) {
            R[i][j] = (i == j) ? 1 : 0;
        }
        t[i] = 0;
    }

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        Point A_transformed[MAX_POINTS];
        for (size_t i = 0; i < sizeA; ++i) {
            matrixVectorMultiply(R, A[i].coords.data(), A_transformed[i].coords.data());
            for (size_t j = 0; j < K; ++j) {
                A_transformed[i].coords[j] += t[j];
            }
        }

        size_t valid_pair_count = 0;
        Point A_valid[MAX_POINTS];
        Point B_valid[MAX_POINTS];
        for (size_t i = 0; i < sizeA; ++i) {
            double best_dist;
            Point nearest_point = tree.nearest(A_transformed[i], best_dist);
            if (best_dist < DISTANCE_THRESHOLD) {
                A_valid[valid_pair_count] = A[i];
                B_valid[valid_pair_count] = nearest_point;
                valid_pair_count++;
            }
        }

        if (valid_pair_count < 3) {
            std::cout << "Not enough valid pairs" << std::endl;
            break;
        }

        double centroid_A[K], centroid_B[K];
        computeCentroid(A_valid, valid_pair_count, centroid_A);
        computeCentroid(B_valid, valid_pair_count, centroid_B);

        subtractCentroid(A_valid, valid_pair_count, centroid_A);
        subtractCentroid(B_valid, valid_pair_count, centroid_B);

        double H[K][K] = {0};
        for (size_t i = 0; i < valid_pair_count; ++i) {
            for (size_t j = 0; j < K; ++j) {
                for (size_t k = 0; k < K; ++k) {
                    H[j][k] += A_valid[i].coords[j] * B_valid[i].coords[k];
                }
            }
        }

        double U[K][K], S[K], V[K][K];
        computeSVD(H, U, S, V);

        double Vt[K][K];
        matrixTranspose(V, Vt);

        double R_new[K][K];
        matrixMultiply(Vt, U, R_new);

        if (R_new[0][0] * (R_new[1][1] * R_new[2][2] - R_new[1][2] * R_new[2][1]) -
            R_new[0][1] * (R_new[1][0] * R_new[2][2] - R_new[1][2] * R_new[2][0]) +
            R_new[0][2] * (R_new[1][0] * R_new[2][1] - R_new[1][1] * R_new[2][0]) < 0) {
            V[2][0] *= -1;
            V[2][1] *= -1;
            V[2][2] *= -1;
            matrixMultiply(Vt, U, R_new);
        }

        double t_new[K];
        double R_centroid_A[K];
        matrixVectorMultiply(R_new, centroid_A, R_centroid_A);
        for (size_t i = 0; i < K; ++i) {
            t_new[i] = centroid_B[i] - R_centroid_A[i];
        }

        double R_diff[K][K];
        for (size_t i = 0; i < K; ++i) {
            for (size_t j = 0; j < K; ++j) {
                R_diff[i][j] = R_new[i][j] - R[i][j];
            }
        }

        double norm_R_diff = 0;
        for (size_t i = 0; i < K; ++i) {
            for (size_t j = 0; j < K; ++j) {
                norm_R_diff += R_diff[i][j] * R_diff[i][j];
            }
        }
        norm_R_diff = std::sqrt(norm_R_diff);

        double norm_t_diff = 0;
        for (size_t i = 0; i < K; ++i) {
            norm_t_diff += (t_new[i] - t[i]) * (t_new[i] - t[i]);
        }
        norm_t_diff = std::sqrt(norm_t_diff);

        if (norm_R_diff < tolerance && norm_t_diff < tolerance) {
            break;
        }

        for (size_t i = 0; i < K; ++i) {
            for (size_t j = 0; j < K; ++j) {
                R[i][j] = R_new[i][j];
            }
            t[i] = t_new[i];
        }
    }
}

int main() {
    Point A[] = {
        { {1, 2, 3} }, { {4, 5, 6} }, { {7, 8, 9} },
        { {2, 3, 4} }, { {5, 6, 7} }, { {8, 9, 10} },
        { {3, 4, 5} }, { {6, 7, 8} }, { {9, 10, 11} }, { {10, 11, 12} }
    };

    Point B[] = {
        { {1.2, 2.1, 2.9} }, { {4.1, 5.0, 6.2} }, { {6.9, 8.2, 9.1} },
        { {2.2, 3.1, 4.3} }, { {5.1, 6.0, 7.2} }, { {8.2, 9.1, 10.0} },
        { {3.1, 4.0, 5.3} }, { {6.0, 7.2, 8.1} }, { {9.2, 10.1, 10.9} },
        { {10.0, 11.1, 12.0} }, { {13.0, 14.1, 15.2} }, { {16.0, 17.1, 18.0} }
    };

    size_t sizeA = sizeof(A) / sizeof(A[0]);
    size_t sizeB = sizeof(B) / sizeof(B[0]);

    double R[K][K];
    double t[K];
    computeICP(A, B, sizeA, sizeB, R, t);

    std::cout << "Rotation matrix R:\n";
    for (size_t i = 0; i < K; ++i) {
        for (size_t j = 0; j < K; ++j) {
            std::cout << R[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Translation vector t:\n";
    for (size_t i = 0; i < K; ++i) {
        std::cout << t[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}

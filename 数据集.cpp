#include <iostream>
#include <cmath>
#include <random>

const int MAX_POINTS = 1000; // 最大点数
const int K = 3; // 空间维度

// 点的类定义
class Point {
public:
    int id;
    double x, y, z; // 坐标
    int cluster; // 聚类标签

    Point(int _id, double _x, double _y, double _z) : id(_id), x(_x), y(_y), z(_z), cluster(-1) {}

    // 计算点之间的空间距离
    double distance(const Point& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        double dz = z - other.z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }
};

// 生成测试数据的类定义
class TestDataGenerator {
public:
    TestDataGenerator(double threshold, int numNoise, double mean, double stddev) : threshold(threshold), numNoise(numNoise), mean(mean), stddev(stddev) {}

    // 生成测试数据集
    void generateData(Point* points[]) {
        generateTargets(points);
        generateNoise(points);
    }

private:
    double threshold; // 目标点之间的最小距离阈值
    int numNoise; // 噪声点的数量
    double mean; // 噪声点坐标的正态分布均值
    double stddev; // 噪声点坐标的正态分布标准差

    // 生成目标点
    void generateTargets(Point* points[]) {
        // 在平面内随机生成 5 个目标点
        for (int i = 0; i < 5; ++i) {
            double x = getRandomNumber(0.0, 10.0);
            double y = getRandomNumber(0.0, 10.0);
            double z = getRandomNumber(0.0, 10.0);
            points[i] = new Point(i + 1, x, y, z);
        }

        // 检查目标点之间的距离，如果小于阈值则重新生成
        for (int i = 0; i < 5; ++i) {
            for (int j = i + 1; j < 5; ++j) {
                if (points[i]->distance(*points[j]) < threshold) {
                    // 重新生成点的坐标
                    double x = getRandomNumber(0.0, 10.0);
                    double y = getRandomNumber(0.0, 10.0);
                    double z = getRandomNumber(0.0, 10.0);
                    points[j]->x = x;
                    points[j]->y = y;
                    points[j]->z = z;
                    // 重新计算距离
                    j = i; // 重新检查目标点之间的距离
                }
            }
        }
    }

    // 生成噪声点
    void generateNoise(Point* points[]) {
        // 在平面内随机生成 numNoise 个噪声点
        for (int i = 0; i < numNoise; ++i) {
            double x = getRandomNormal(mean, stddev);
            double y = getRandomNormal(mean, stddev);
            double z = getRandomNormal(mean, stddev);
            points[i + 5] = new Point(-1, x, y, z); // 噪声点的 ID 设为 -1
        }
    }

    // 生成指定范围内的随机数
    double getRandomNumber(double min, double max) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(min, max);
        return dis(gen);
    }

    // 生成正态分布的随机数
    double getRandomNormal(double mean, double stddev) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<> dis(mean, stddev);
        return dis(gen);
    }
};

int main() {
    Point* points[MAX_POINTS];

    // 创建测试数据生成器，设置阈值、噪声点数量、均值和标准差
    TestDataGenerator dataGenerator(2.0, 10, 5.0, 1.0);
    // 生成测试数据集
    dataGenerator.generateData(points);

    // 输出生成的测试数据
    for (int i = 0; i < MAX_POINTS; ++i) {
        if (points[i] != nullptr) {
            std::cout << "Point " << points[i]->id << ": (" << points[i]->x << ", " << points[i]->y << ", " << points[i]->z << ")" << std::endl;
        }
    }

    // 释放动态分配的内存
    for (int i = 0; i < MAX_POINTS; ++i) {
        delete points[i];
    }

    return 0;
}
//静态
#include <iostream>
#include <cmath>
#include <random>

const int MAX_POINTS = 1000; // 最大点数
const int K = 3; // 空间维度

// 点的类定义
class Point {
public:
    int id;
    double x, y, z; // 坐标
    int cluster; // 聚类标签

    Point(int _id, double _x, double _y, double _z) : id(_id), x(_x), y(_y), z(_z), cluster(-1) {}

    // 计算点之间的空间距离
    double distance(const Point& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        double dz = z - other.z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }
};

// 生成测试数据的类定义
class TestDataGenerator {
public:
    TestDataGenerator(double threshold, int numNoise, double mean, double stddev) : threshold(threshold), numNoise(numNoise), mean(mean), stddev(stddev) {}

    // 生成测试数据集
    void generateData(Point points[]) {
        generateTargets(points);
        generateNoise(points);
    }

private:
    double threshold; // 目标点之间的最小距离阈值
    int numNoise; // 噪声点的数量
    double mean; // 噪声点坐标的正态分布均值
    double stddev; // 噪声点坐标的正态分布标准差

    // 生成目标点
    void generateTargets(Point points[]) {
        // 在平面内随机生成 5 个目标点
        for (int i = 0; i < 5; ++i) {
            double x = getRandomNumber(0.0, 10.0);
            double y = getRandomNumber(0.0, 10.0);
            double z = getRandomNumber(0.0, 10.0);
            points[i] = Point(i + 1, x, y, z);
        }

        // 检查目标点之间的距离，如果小于阈值则重新生成
        for (int i = 0; i < 5; ++i) {
            for (int j = i + 1; j < 5; ++j) {
                if (points[i].distance(points[j]) < threshold) {
                    // 重新生成点的坐标
                    double x = getRandomNumber(0.0, 10.0);
                    double y = getRandomNumber(0.0, 10.0);
                    double z = getRandomNumber(0.0, 10.0);
                    points[j] = Point(j + 1, x, y, z);
                    // 重新计算距离
                    j = i; // 重新检查目标点之间的距离
                }
            }
        }
    }

    // 生成噪声点
    void generateNoise(Point points[]) {
        // 在平面内随机生成 numNoise 个噪声点
        for (int i = 0; i < numNoise; ++i) {
            double x = getRandomNormal(mean, stddev);
            double y = getRandomNormal(mean, stddev);
            double z = getRandomNormal(mean, stddev);
            points[i + 5] = Point(-1, x, y, z); // 噪声点的 ID 设为 -1
        }
    }

    // 生成指定范围内的随机数
    double getRandomNumber(double min, double max) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(min, max);
        return dis(gen);
    }

    // 生成正态分布的随机数
    double getRandomNormal(double mean, double stddev) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<> dis(mean, stddev);
        return dis(gen);
    }
};

int main() {
    Point points[MAX_POINTS];

    // 创建测试数据生成器，设置阈值、噪声点数量、均值和标准差
    TestDataGenerator dataGenerator(2.0, 10, 5.0, 1.0);
    // 生成测试数据集
    dataGenerator.generateData(points);

    // 输出生成的测试数据
    for (int i = 0; i < MAX_POINTS; ++i) {
        if (points[i].id != 0) {
            std::cout << "Point " << points[i].id << ": (" << points[i].x << ", " << points[i].y << ", " << points[i].z << ")" << std::endl;
        }
    }

    return 0;
}

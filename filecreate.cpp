#include <iostream>
#include <string>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#endif

bool createDirectory(const std::string& path) {
#ifdef _WIN32
    // 在Windows中创建目录
    if (CreateDirectory(path.c_str(), nullptr) || GetLastError() == ERROR_ALREADY_EXISTS) {
        return true; // 目录创建成功或已存在
    } else {
        std::cerr << "Error creating directory: " << path << " (Error code: " << GetLastError() << ")\n";
        return false; // 创建目录失败
    }
#else
    // 在Linux中创建目录
    if (mkdir(path.c_str(), 0755) == 0 || errno == EEXIST) {
        return true; // 目录创建成功或已存在
    } else {
        std::cerr << "Error creating directory: " << path << " (Error code: " << errno << ")\n";
        return false; // 创建目录失败
    }
#endif
}

int main() {
    std::string path = "d:/Event-driven"; // 你想要创建的路径
    if (createDirectory(path)) {
        std::cout << "Directory created successfully or already exists: " << path << std::endl;
    } else {
        std::cout << "Failed to create directory: " << path << std::endl;
    }
    return 0;
}

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

struct Point3 {
    double x, y, z;
};

int main() {
    std::vector<Point3> vertices_mesh;
    std::vector<std::vector<int>> faces_mesh;

    std::ifstream obj_file("BaseModels.obj");
    std::string line;
    while (std::getline(obj_file, line)) {
        std::istringstream iss(line);
        char type;
        iss >> type;
        if (type == 'v') {
            // Parse vertex coordinates
            Point3 vertex;
            iss >> vertex.x >> vertex.y >> vertex.z;
            vertices_mesh.push_back(vertex);
        } else if (type == 'f') {
            // Parse face vertex indices
            std::vector<int> face;
            int index;
            while (iss >> index) {
                face.push_back(index - 1); // OBJ indices are 1-based
                iss.ignore(256, ' '); // ignore texture and normal indices
            }
            faces_mesh.push_back(face);
        }
    }

    // Output vertices
    std::cout << "auto vertices_mesh = std::vector<Point3>{";
    for (const auto& vertex : vertices_mesh) {
        std::cout << "{" << vertex.x << ", " << vertex.y << ", " << vertex.z << "}, ";
    }
    std::cout << "};\n";

    // Output faces
    std::cout << "auto faces_mesh = std::vector<std::vector<int>>{";
    for (const auto& face : faces_mesh) {
        std::cout << "{";
        for (const auto& index : face) {
            std::cout << index << ", ";
        }
        std::cout << "}, ";
    }
    std::cout << "};\n";

    return 0;
}

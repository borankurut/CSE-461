#include "parser/pugixml/src/pugixml.hpp"
#include <iostream>
#include <vector>
#include <sstream>
#include <tuple>

struct Vector3 {
    float x, y, z;
};

struct Camera {
    Vector3 position;
    Vector3 gaze;
    Vector3 up;
    std::tuple<float, float, float, float> nearPlane;
    float nearDistance;
    std::tuple<int, int> imageResolution;
};

struct Light {
    std::string type;
    int id;
    Vector3 position;
    Vector3 intensity;
};

struct Material {
    int id;
    Vector3 ambient;
    Vector3 diffuse;
    Vector3 specular;
    Vector3 mirrorReflectance;
    float phongExponent;
};

struct Face {
    std::vector<int> vertices;
};

struct Mesh {
    int id;
    int materialId;
    std::vector<Face> faces;
};

int main() {
    pugi::xml_document doc;
    if (!doc.load_file("scene.xml")) {
        std::cerr << "Error loading XML file." << std::endl;
        return 1;
    }

    // Extract camera data
    Camera camera;
    auto camera_node = doc.child("scene").child("camera");
    camera.position = { camera_node.child("position").text().as_float(),
                        camera_node.child("position").next_sibling("gaze").text().as_float(),
                        camera_node.child("position").next_sibling("gaze").next_sibling("up").text().as_float() };
    camera.gaze = { camera_node.child("gaze").text().as_float(),
                    camera_node.child("gaze").next_sibling("up").text().as_float(),
                    camera_node.child("gaze").next_sibling("up").next_sibling("nearPlane").text().as_float() };
    camera.up = { camera_node.child("up").text().as_float(),
                  camera_node.child("up").next_sibling("nearPlane").text().as_float(),
                  camera_node.child("up").next_sibling("nearPlane").next_sibling("neardistance").text().as_float() };
    std::istringstream near_plane_ss(camera_node.child("nearPlane").text().as_string());
    near_plane_ss >> std::get<0>(camera.nearPlane) >> std::get<1>(camera.nearPlane) >> std::get<2>(camera.nearPlane) >> std::get<3>(camera.nearPlane);
    camera.nearDistance = camera_node.child("neardistance").text().as_float();
    std::istringstream image_resolution_ss(camera_node.child("imageresolution").text().as_string());
    image_resolution_ss >> std::get<0>(camera.imageResolution) >> std::get<1>(camera.imageResolution);

    // Extract lights
    std::vector<Light> lights;
    for (pugi::xml_node light_node = doc.child("scene").child("lights").first_child(); light_node; light_node = light_node.next_sibling()) {
        Light light;
        light.type = light_node.name();
        if (light.type == "ambientlight") {
            light.intensity = { light_node.text().as_float(), 0.0f, 0.0f };
        }
        else {
            light.id = light_node.attribute("id").as_int();
            light.position = { light_node.child("position").text().as_float(),
                                light_node.child("position").next_sibling("intensity").text().as_float(),
                                light_node.child("position").next_sibling("intensity").next_sibling().text().as_float() };
            light.intensity = { light_node.child("intensity").text().as_float(),
                                light_node.child("intensity").next_sibling().text().as_float(),
                                light_node.child("intensity").next_sibling().next_sibling().text().as_float() };
        }
        lights.push_back(light);
    }

    // Extract materials
    std::vector<Material> materials;
    for (pugi::xml_node material_node = doc.child("scene").child("materials").first_child(); material_node; material_node = material_node.next_sibling()) {
        Material material;
        material.id = material_node.attribute("id").as_int();
        material.ambient = { material_node.child("ambient").text().as_float(),
                            material_node.child("ambient").next_sibling("diffuse").text().as_float(),
                            material_node.child("ambient").next_sibling("diffuse").next_sibling("specular").text().as_float() };
        material.diffuse = { material_node.child("diffuse").text().as_float(),
                            material_node.child("diffuse").next_sibling("specular").text().as_float(),
                            material_node.child("diffuse").next_sibling("specular").next_sibling("mirrorreflactance").text().as_float() };
        material.specular = { material_node.child("specular").text().as_float(),
                             material_node.child("specular").next_sibling("mirrorreflactance").text().as_float(),
                             material_node.child("specular").next_sibling("mirrorreflactance").next_sibling("phongexponent").text().as_float() };
        material.mirrorReflectance = { material_node.child("mirrorreflactance").text().as_float(),
                                      material_node.child("mirrorreflactance").next_sibling("phongexponent").text().as_float(),
                                      0.0f };
        material.phongExponent = material_node.child("phongexponent").text().as_float();
        materials.push_back(material);
    }

    // Extract vertex data
    std::vector<Vector3> vertices;
    std::istringstream vertex_data_ss(doc.child("scene").child("vertexdata").child("data").text().as_string());
    while (!vertex_data_ss.eof()) {
        Vector3 vertex;
        vertex_data_ss >> vertex.x >> vertex.y >> vertex.z;
        vertices.push_back(vertex);
    }

    // Extract objects
    std::vector<Mesh> meshes;
    for (pugi::xml_node mesh_node = doc.child("scene").child("objects").first_child(); mesh_node; mesh_node = mesh_node.next_sibling()) {
        Mesh mesh;
        mesh.id = mesh_node.attribute("id").as_int();
        mesh.materialId = mesh_node.child("materialid").text().as_int();
        for (pugi::xml_node face_node = mesh_node.child("faces").first_child(); face_node; face_node = face_node.next_sibling()) {
            Face face;
            std::istringstream face_ss(face_node.text().as_string());
            while (!face_ss.eof()) {
                int vertex_index;
                face_ss >> vertex_index;
                face.vertices.push_back(vertex_index);
            }
            mesh.faces.push_back(face);
        }
        meshes.push_back(mesh);
    }

    // Output extracted data
    std::cout << "Camera Data:" << std::endl;
    std::cout << "  Position: " << camera.position.x << " " << camera.position.y << " " << camera.position.z << std::endl;
    std::cout << "  Gaze: " << camera.gaze.x << " " << camera.gaze.y << " " << camera.gaze.z << std::endl;
    std::cout << "  Up: " << camera.up.x << " " << camera.up.y << " " << camera.up.z << std::endl;
    std::cout << "  Near Plane: " << std::get<0>(camera.nearPlane) << " " << std::get<1>(camera.nearPlane) << " " << std::get<2>(camera.nearPlane) << " " << std::get<3>(camera.nearPlane) << std::endl;
    std::cout << "  Near Distance: " << camera.nearDistance << std::endl;
    std::cout << "  Image Resolution: " << std::get<0>(camera.imageResolution) << "x" << std::get<1>(camera.imageResolution) << std::endl;

    std::cout << "\nLight Sources:" << std::endl;
    for (const auto& light : lights) {
        std::cout << "  Type: " << light.type << std::endl;
        if (light.type != "ambientlight") {
            std::cout << "  ID: " << light.id << std::endl;
            std::cout << "  Position: " << light.position.x << " " << light.position.y << " " << light.position.z << std::endl;
        }
        std::cout << "  Intensity: " << light.intensity.x << " " << light.intensity.y << " " << light.intensity.z << std::endl;
    }

    std::cout << "\nMaterials:" << std::endl;
    for (const auto& material : materials) {
        std::cout << "  ID: " << material.id << std::endl;
        std::cout << "  Ambient: " << material.ambient.x << " " << material.ambient.y << " " << material.ambient.z << std::endl;
        std::cout << "  Diffuse: " << material.diffuse.x << " " << material.diffuse.y << " " << material.diffuse.z << std::endl;
        std::cout << "  Specular: " << material.specular.x << " " << material.specular.y << " " << material.specular.z << std::endl;
        std::cout << "  Mirror Reflectance: " << material.mirrorReflectance.x << " " << material.mirrorReflectance.y << " " << material.mirrorReflectance.z << std::endl;
        std::cout << "  Phong Exponent: " << material.phongExponent << std::endl;
    }

    std::cout << "\nVertex Data:" << std::endl;
    for (const auto& vertex : vertices) {
        std::cout << "  " << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
    }

    std::cout << "\nObjects:" << std::endl;
    for (const auto& mesh : meshes) {
        std::cout << "  Mesh ID: " << mesh.id << std::endl;
        std::cout << "  Material ID: " << mesh.materialId << std::endl;
        std::cout << "  Faces:" << std::endl;
        for (const auto& face : mesh.faces) {
            std::cout << "    ";
            for (int vertex_index : face.vertices) {
                std::cout << vertex_index << " ";
            }
            std::cout << std::endl;
        }
    }

    return 0;
}

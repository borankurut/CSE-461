#include "parser/pugixml/src/pugixml.hpp"
#include "color.h"
#include "camera.h"
#include "directional_light.h"
#include "material.h"
#include "mesh_object.h"
#include "point_light.h"
#include "raytracing.h"
#include "light.h"
#include "hittable.h"
#include "hittable_list.h"
#include "sphere.h"
#include "vec3.h"
#include <iostream>
#include <memory>
#include "ambient_light.h"
#include <ostream>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <sstream>
#include <chrono>

std::mutex mtx; // Mutex for synchronization

struct World{
	Hittable_list objects;
	Color background;
    std::vector<shared_ptr<Light>> lights;
	shared_ptr<Light> ambient_light;
	int max_ray_depth;
};

void addMatterial(Hittable & hittable, std::string id, const std::vector<shared_ptr<Material>> &material_list){
	for(const auto & material: material_list){
		if(id == material->id){
			hittable.material = material;
			return;
		}
	}

	std::clog << "ERROR: Material list doesn't have any material: " << id << "\n";
}

Color shadeRay(const Ray& r, int depth, const World& world){
	if(depth <= 0)
		return Color(0,0,0);

    HitRecord rec;

	bool is_hit = world.objects.hit(r, 0, infinity, rec);

    if (is_hit) {
		Color c = world.ambient_light->illuminate(r, rec);

		for (auto light : world.lights){
			Vec3 to_light = (light->dir(rec.p)); 
			//std::clog << to_light;
			auto distance = to_light.length();
			to_light = unit_vector(to_light);

			auto epsilon = 0.0000001;

			HitRecord shadow_rec;

			bool is_shadow = world.objects.hit(Ray(rec.p, to_light), epsilon, distance, shadow_rec);

			if(!is_shadow){
				c += light->illuminate(r, rec);
				fix_color(c);

				Color km = rec.hit_surface_p->material->km;

				if(!(km.x() == 0 && km.y() == 0 && km.z() == 0)){
					Vec3 reflectRay_dir = reflect(-to_light, rec.normal);
					Ray reflectRay = Ray(rec.p + epsilon * reflectRay_dir, reflectRay_dir);
					Color reflectColor = shadeRay(reflectRay, depth - 1, world) * km;
					c = 0.75 * c + 0.25 * reflectColor;
				}
			}
		}

		return c;
    }
	return world.background;
}

Color shadeRayShowLight(Point3 lp, const Ray& r, const World& world){ // debug purpose
	Sphere s = Sphere(lp, 0.01);
	HitRecord rec; 
	bool is_hit_light = s.hit(r, 0, infinity, rec);
	if(is_hit_light)
		return Color(1, 1, 1);
	else 
		return shadeRay(r, 10, world);

}

void processChunk(int start, int end, int ny, int nx, const Camera& camera, const World& world, std::vector<std::vector<Color>>& pixel_colors) {
    for (int j = start; j < end; j++) {
        for (int i = nx - 1; i > 0; i--) {
            const auto& np = camera.nearplane();
            auto u = np.left + (np.right - np.left) * (i + 0.5) / nx;
            auto v = np.bottom + (np.top - np.bottom) * (j + 0.5) / ny;
            Vec3 ray_o = camera.position();
            Vec3 ray_d = (u * camera.u()) + (v * camera.v()) - camera.w() * camera.neardistance();
            Ray r = Ray(ray_o, ray_d);
            Color pixel_color = shadeRay(r,world.max_ray_depth, world);

            // Lock the mutex before writing to std::cout
            std::lock_guard<std::mutex> lock(mtx);
			pixel_colors[j][i] = pixel_color;
//            write_color(std::cout, pixel_color);
        }
    }
	std::clog<< "Chunk from " << start << " to " << end << " is finished.\n";

}

Vec3 stringToVector(const std::string& str){
	std::istringstream streamm(str);

	double x, y, z; 
	streamm >> x>> y >> z;
	return Vec3(x, y, z);
}


shared_ptr<PointLight> parsePointLight(const pugi::xml_node& lightNode);
shared_ptr<DirectionalLight> parseTriLight(const pugi::xml_node& lightNode);
shared_ptr<Material> parseMaterial(const pugi::xml_node& materialNode);
void parseVertex(const pugi::xml_node& materialNode, std::vector<Point3>& vertices);
shared_ptr<MeshObject> parseMesh(const pugi::xml_node& objectNode, const std::vector<Point3>& vertices);

int main() {
    // World
	World world;
	std::string filename = "scene_what.xml";

    pugi::xml_document doc;
    if (!doc.load_file(filename.c_str())) {
        std::cerr << "Failed to load XML file: " << filename << std::endl;
        return 0;
    }

    // Extract maxRayDepth from XML
    world.max_ray_depth = doc.child("scene").child("maxraytracedepth").text().as_int();

    // Extract background color from XML
    // Assuming background color format is "R G B"
	std::string background_str = doc.child("scene").child("backgroundColor").text().as_string();

	world.background = stringToVector(background_str);
	world.background.set_x(world.background.x()/255);world.background.set_y(world.background.y()/255);world.background.set_z(world.background.z()/255);

	//Camera
	Vec3 camera_pos = stringToVector(doc.child("scene").child("camera").child("position").text().as_string());
	Vec3 camera_gaze = stringToVector(doc.child("scene").child("camera").child("gaze").text().as_string());
	Vec3 camera_up = stringToVector(doc.child("scene").child("camera").child("up").text().as_string());

	std::string camera_nearplane = doc.child("scene").child("camera").child("nearPlane").text().as_string();
	std::istringstream nps(camera_nearplane);

	double camera_near_l, camera_near_r, camera_near_b, camera_near_t;
	nps >> camera_near_l >> camera_near_r >> camera_near_b >> camera_near_t;
	double neardistance = doc.child("scene").child("camera").child("neardistance").text().as_double();

	std::istringstream camera_res(doc.child("scene").child("camera").child("imageresolution").text().as_string());
	int res_x, res_y;
	camera_res >> res_x >> res_y;


	Camera camera = Camera(camera_pos, camera_gaze, camera_up, Nearplane(camera_near_l, camera_near_r, camera_near_b, camera_near_t), neardistance, res_x, res_y);
	/* std::clog << camera_pos << std::endl; */
	/* std::clog << camera_gaze<< std::endl; */
	/* std::clog << camera_up<< std::endl; */
	/* std::clog << camera_near_l<< std::endl; */
	/* std::clog << camera_near_r<< std::endl; */
	/* std::clog << camera_near_b<< std::endl; */
	/* std::clog << camera_near_t<< std::endl; */
	/* std::clog << neardistance<< std::endl; */
	/* std::clog << res_x<< std::endl; */
	/* std::clog << res_y<< std::endl; */

	Vec3 ambvec= stringToVector(doc.child("scene").child("lights").child("ambientlight").text().as_string());
	world.ambient_light = make_shared<AmbientLight>(AmbientLight::from_rgb(ambvec));

	for(auto lightNode : doc.child("scene").child("lights").children()){
		std::string name = lightNode.name();
		if(name == "pointlight"){
			auto pl = parsePointLight(lightNode);
			world.lights.push_back(pl);
		}
		else if(name == "triangularlight"){
			auto tl = parseTriLight(lightNode);
			world.lights.push_back(tl);
		}
	}

	std::vector<shared_ptr<Material>> material_list;

	for(auto material_node: doc.child("scene").child("materials").children()){
		material_list.push_back(parseMaterial(material_node));
	}

	std::vector<Point3> vertices;
	parseVertex(doc.child("scene").child("vertexdata"), vertices);

	for(auto objectNode: doc.child("scene").child("objects").children()){
		std::string material_id = objectNode.child("materialid").text().as_string();
		auto to_add = parseMesh(objectNode, vertices);

		addMatterial(*to_add, material_id, material_list);

		world.objects.add(to_add);
	}


 // Render
	int nx = camera.nx();
	int ny = camera.ny();
    std::cout << "P3\n" << nx << ' ' << ny << "\n255\n";

	std::vector<std::vector<Color>> pixel_colors(ny+1, std::vector<Color>(nx+1));
	auto start_rt_time = std::chrono::system_clock::now();

	unsigned int num_threads = std::thread::hardware_concurrency();
	std::clog << "num_threads: " << num_threads << "\n";
	int chunk_size = ny / num_threads;

	// create threads by dividing ny to thread amount.
	std::vector<std::thread> threads;
	for (unsigned int t = 0; t < num_threads; ++t) {
		int start = t * chunk_size;
		int end = (t == num_threads - 1) ? ny : (t + 1) * chunk_size;
		threads.emplace_back(processChunk, start, end, ny, nx, std::ref(camera), std::ref(world), std::ref(pixel_colors));
	}

	// join threads
	for (auto& thread : threads) {
		thread.join();
	}


	auto end_rt_time = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end_rt_time - start_rt_time;
	std::clog << "Ray tracing finished in " << elapsed_seconds.count() << " s." << "\n";

	std::clog << "Writing into ppm file...\n";

	for (int j = 0; j < ny; ++j) {
        for (int i = nx; i > 0; --i) {
            mtx.lock(); // Lock mutex before accessing shared resource
            write_color(std::cout, pixel_colors[j][i]); // Output color value to std::cout
            mtx.unlock(); // Unlock mutex after accessing shared resource
        }
        std::cout << "\n"; // Newline after each row of pixels
    }

	std::clog << "Writing finished.\n";





	return 0;
}

shared_ptr<PointLight> parsePointLight(const pugi::xml_node& lightNode){
	auto pos = stringToVector(lightNode.child("position").text().as_string());
	auto intensity = stringToVector(lightNode.child("intensity").text().as_string());

	return 	make_shared<PointLight>(PointLight::fromIntensity(pos, intensity));

}

shared_ptr<DirectionalLight> parseTriLight(const pugi::xml_node& lightNode){
	auto v1 = stringToVector(lightNode.child("vertex1").text().as_string());
	auto v2 = stringToVector(lightNode.child("vertex2").text().as_string());
	auto v3 = stringToVector(lightNode.child("vertex3").text().as_string());

	auto intensity = stringToVector(lightNode.child("intensity").text().as_string());

	return make_shared<DirectionalLight>(DirectionalLight::fromTriangle(v1, v2, v3, intensity));
}
shared_ptr<Material> parseMaterial(const pugi::xml_node& materialNode){
	std::string id = materialNode.attribute("id").as_string();

	auto ka = stringToVector(materialNode.child("ambient").text().as_string());
	auto kd = stringToVector(materialNode.child("diffuse").text().as_string());
	auto ks = stringToVector(materialNode.child("specular").text().as_string());
	auto km = stringToVector(materialNode.child("mirrorreflactance").text().as_string());
	double ph = materialNode.child("phongexponent").text().as_double();

	return make_shared<Material>(Material(id, ka, kd, ks, km, ph));
}


void parseVertex(const pugi::xml_node& vertexNode, std::vector<Point3>& vertices){
	std::string str = vertexNode.text().as_string();
	std::istringstream streamm(str);

	double x, y, z; 
	while(streamm >> x >> y >> z){
		vertices.push_back(Vec3(x, y, z));
	}
}


shared_ptr<MeshObject> parseMesh(const pugi::xml_node& objectNode, const std::vector<Point3>& vertices){
	std::vector<std::vector<int>> faces;

	std::string str = objectNode.child("faces").text().as_string();
	std::istringstream streamm(str);

	int x, y, z; 
	while(streamm >> x >> y >> z){
		faces.push_back(std::vector<int>({x, y, z}));
	}

	return make_shared<MeshObject>(MeshObject(vertices, faces));
}


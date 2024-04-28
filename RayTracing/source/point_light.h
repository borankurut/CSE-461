#ifndef POINT_LIGHT_H
#define POINT_LIGHT_H

#include "color.h"
#include "light.h"
#include "vec3.h"
#include <cmath>
#include <iostream>

class Point_light : public Light{
  
public:
	Point_light(const Color& c, const Point3& p) : c(c), p(p) {}

	// I assumed the intensity values are between 0 and 1000.
	static Point_light fromIntensity(const Point3& p, const Vec3& intensity){

		return Point_light(Color(intensity.x()/1000.0, intensity.y()/1000.0, intensity.z()/1000.0), p);

	}

	Vec3 dir(const Point3& point) override {
		return p - point;
	}

	Color illuminate(const Ray& r, const Hit_record& rec) override{
		// Point3 hit_point = rec.p;
		// Vec3 p_min_x = p - hit_point;
		// auto p_min_x_len = p_min_x.length();
		//
		// // Check if the distance is zero to avoid division by zero
		// if (p_min_x_len == 0.0) {
		// 	return Color(1, 0, 0); // Or handle this case appropriately
		// }
		//
		// Vec3 light_direction = p_min_x / p_min_x_len;
		// Vec3 n = rec.normal;
		//
		// // Compute diffuse lighting contribution
		// Color diffuse = fmax(0.0, dot(n, light_direction)) * c / (p_min_x_len * p_min_x_len);
		//
		// // Access material if it exists
		// if (rec.hit_surface_p && rec.hit_surface_p->material) {
		// 	Vec3 v = -r.direction();
		// 	Color specular = rec.hit_surface_p->material->evaluate(light_direction, v, n);
		// 	return diffuse * specular;
		// } 
		// else 
		// 	return diffuse;
		
		Color E = Color(0, 0, 0);
		Color k = Color(0, 0, 0);

		Point3 hit_point = rec.p;

		Vec3 p_min_x = p - hit_point;

		auto p_min_x_len = p_min_x.length();

		Vec3 light_direction = p_min_x / p_min_x_len;

		Vec3 n = rec.normal;

		if(p_min_x_len != 0)
			E = std::max(0.0, dot(n, light_direction)) * c / pow(p_min_x_len, 2);

		Vec3 v = -r.direction();

		if(rec.hit_surface_p && rec.hit_surface_p->material)
			k = rec.hit_surface_p->material->evaluate(light_direction, v, n);

		return E*k;
	}
private:
	Color c;
	Point3 p;
};

#endif


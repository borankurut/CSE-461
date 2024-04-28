#ifndef DIRECTIONAL_LIGHT_H
#define DIRECTIONAL_LIGHT_H

#include "light.h"
#include "vec3.h"

class Directional_light: public Light{
  
public:
	Directional_light (const Color& c, const Vec3& d) : c(c), d(d) {}

	/**
	 *<vertex1> x y z </vertex1><vertex2> x y z </vertex2>1
	<vertex3> x y z </vertex3><intensity> x y z </intensity>
	(vertex1-vertex2) x (vertex1-vertex3)

	I assumed the intensity values are between 0 and 1000.
	*/
	static Directional_light fromTriangle(Vec3 vertex1, Vec3 vertex2,Vec3 vertex3, Vec3 intensity){

		Vec3 direction = cross(vertex1 - vertex2, vertex1 - vertex3);
//		intensity.set_x(intensity.x() / 1000.0);
		Vec3 intensityFixed = Vec3(intensity.x()/ 1000.0, intensity.y()/1000.0, intensity.z()/1000.0);

		return Directional_light(intensityFixed, direction);
	}

	Vec3 dir(const Point3& point) override {
		return d;
	}

	Color illuminate(const Ray& r, const Hit_record& rec) override{
		Color E = Color(0, 0, 0);
		Color k = Color(0, 0, 0);

		Vec3 n = rec.normal;

		E = std::max(0.0, dot(n, d)) * c;
		if(dot(n,d) < 0)
			return (Color(0, 0, 0));

		Vec3 v = -r.direction();

		if(rec.hit_surface_p && rec.hit_surface_p->material)
			k = rec.hit_surface_p->material->evaluate(d, v, n);

		return E*k;
	}
private:
	Color c;
	Vec3 d;
};

#endif


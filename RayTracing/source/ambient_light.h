#include "light.h"
#include "vec3.h"
#include <iostream>

class AmbientLight: public Light{
  
public:
	AmbientLight(const Color& c) : c(c){}

	static AmbientLight from_rgb(const Vec3& rgb){
		return AmbientLight(Color(rgb.x()/255.0, rgb.y()/255.0, rgb.z()/255.0));
	}

	Vec3 dir(const Point3& point) override {
		return Vec3(0, 0, 0);
	}

	Color illuminate(const Ray& r, const HitRecord& rec) override{
		Color ka = rec.hit_surface_p->material->ka;
		return c * ka;
	}
private:
	Color c;
};


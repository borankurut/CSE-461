#include "light.h"
#include "vec3.h"
#include <iostream>

class Ambient_light: public Light{
  
public:
	Ambient_light(const Color& c) : c(c){}

	static Ambient_light from_rgb(const Vec3& rgb){
		return Ambient_light(Color(rgb.x()/255.0, rgb.y()/255.0, rgb.z()/255.0));
	}

	Vec3 dir(const Point3& point) override {
		return Vec3(0, 0, 0);
	}

	Color illuminate(const Ray& r, const Hit_record& rec) override{
		Color ka = rec.hit_surface_p->material->ka;
		return c * ka;
	}
private:
	Color c;
};


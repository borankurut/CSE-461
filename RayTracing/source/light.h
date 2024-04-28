#ifndef LIGHT_H
#define LIGHT_H

#include "color.h"
#include "ray.h"
#include "hittable.h"

class Light{
public:
	virtual Color illuminate(const Ray& r, const Hit_record& rec) = 0;
	virtual Vec3 dir(const Point3& point) = 0;
};

#endif


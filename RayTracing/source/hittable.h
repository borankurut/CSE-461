#ifndef HITTABLE_H
#define HITTABLE_H

#include "ray.h"
#include "material.h"
#include <memory>

class Hittable;

class Hit_record {
  public:
    Point3 p;
    Vec3 normal;
    double t;

	const Hittable* hit_surface_p;

    bool front_face;

    void set_face_normal(const Ray& r, const Vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

class Hittable {
  public:
    virtual ~Hittable() = default;

    virtual bool hit(const Ray& r, double ray_tmin, double ray_tmax, Hit_record& rec) const = 0;

	std::shared_ptr<Material> material;
};

#endif

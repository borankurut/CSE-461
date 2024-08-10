#ifndef BOX_H
#define BOX_H

#include "hittable.h"
#include "vec3.h"

class HitBox: public Hittable {
  public:
    HitBox() : _min_point(Point3(0,0,0)), _max_point(Point3(0, 0, 0)) {}

    HitBox(const Point3& min_point, const Point3& max_point) : _min_point(min_point), _max_point(max_point) {}

	// only returns whether the ray hit the box, doesn't fill the record.
    bool hit(const Ray& r, double ray_tmin, double ray_tmax, HitRecord& rec) const override {
        for (int i = 0; i < 3; ++i) {
            double inverted_direction = 1.0 / r.direction()[i];
            double t0 = (_min_point[i] - r.origin()[i]) * inverted_direction;
            double t1 = (_max_point[i] - r.origin()[i]) * inverted_direction;

            if (inverted_direction < 0.0) 
				std::swap(t0, t1);

            ray_tmin = t0 > ray_tmin ? t0 : ray_tmin;
            ray_tmax = t1 < ray_tmax ? t1 : ray_tmax;

            if (ray_tmax <= ray_tmin)
                return false;
        }

        return true;
    }

  private:
    Point3 _min_point;
    Point3 _max_point;
};

#endif


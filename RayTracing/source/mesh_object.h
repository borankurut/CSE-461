#ifndef MESH_OBJECT_H
#define MESH_OBJECT_H

#include "hittable.h"
#include "vec3.h"
#include <vector>

class MeshObject : public Hittable {
public:
    MeshObject(const std::vector<Point3>& vertices, const std::vector<std::vector<int>>& faces) 
        : vertices(vertices), faces(faces) {}

    bool hit(const Ray& r, double ray_tmin, double ray_tmax, Hit_record& rec) const override {
        bool hit_anything = false;
        double closest_so_far = ray_tmax;

        // Check each face
        for (const auto& face : faces) {
            const Point3& v0 = vertices[face[0]];
            const Point3& v1 = vertices[face[1]];
            const Point3& v2 = vertices[face[2]];

            Vec3 e1 = v1 - v0;
            Vec3 e2 = v2 - v0;
            Vec3 h = cross(r.direction(), e2);
            double a = dot(e1, h);

            if (a > -0.000001 && a < 0.000001)
                continue;

            double f = 1 / a;
            Vec3 s = r.origin() - v0;
            double u = f * dot(s, h);

            if (u < 0 || u > 1)
                continue;

            Vec3 q = cross(s, e1);
            double v = f * dot(r.direction(), q);

            if (v < 0 || u + v > 1)
                continue;

            double t = f * dot(e2, q);

            if (t < ray_tmin || t > closest_so_far)
                continue;

            // Hit
            closest_so_far = t;
            hit_anything = true;

            // Record hit information
            rec.t = t;
            rec.p = r.at(t);
            Vec3 outward_normal = unit_vector(cross(e1, e2)); // Assuming counterclockwise winding order
            rec.set_face_normal(r, outward_normal);
            rec.hit_surface_p = this;
        }

        return hit_anything;
    }

private:
    std::vector<Point3> vertices;
    std::vector<std::vector<int>> faces;
};

#endif


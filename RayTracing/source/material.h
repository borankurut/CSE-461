#ifndef MATERIAL_H 
#define MATERIAL_H 

#include "color.h"
#include "vec3.h"
#include <string>

class Material{
public:
	Material(std::string id, Color ka, Color kd, Color ks, Color km, double phongeE): id(id), ka(ka), kd(kd), ks(ks), km(km), phongeE(phongeE){}

	std::string id;
	Color ka;
	Color kd;
	Color ks;
	Color km;

	double phongeE;

	Color evaluate(Vec3 l, Vec3 v, Vec3 n){
		double diffuse_intensity = std::max(0.0, dot(n, l));

		Color diffuse_term = kd * diffuse_intensity;

		Vec3 h = unit_vector(l + v);

		double specular_intensity = std::pow(std::max(0.0, dot(n, h)), phongeE);

		Color specular_term = ks * specular_intensity;

		// Combine ambient, diffuse, and specular terms
		Color total_reflectance = (1.0/3.0) * ka + (1.0/3.0) * diffuse_term + (1.0/3.0) *  specular_term;

		return total_reflectance;
	}

};

#endif

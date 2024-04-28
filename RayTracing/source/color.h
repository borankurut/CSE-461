#ifndef COLOR_H
#define COLOR_H

#include "vec3.h"

#include <iostream>

using Color = Vec3;

void write_color(std::ostream& out, const Color& pixel_color) {
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // Translate the [0,1] component values to the byte range [0,255].
    int rbyte = int(255.999 * r);
    int gbyte = int(255.999 * g);
    int bbyte = int(255.999 * b);

    // Write out the pixel color components.
    out << rbyte << ' ' << gbyte << ' ' << bbyte << '\n';
}

void fix_color(Color& pixel_color){
	if(pixel_color.x() > 1)
		pixel_color.set_x(1);
	else if(pixel_color.x() < 0)
		pixel_color.set_x(0);

	if(pixel_color.y() > 1)
		pixel_color.set_y(1);
	else if(pixel_color.y() < 0)
		pixel_color.set_y(0);

	if(pixel_color.z() > 1)
		pixel_color.set_z(1);
	else if(pixel_color.z() < 0)
		pixel_color.set_z(0);
}

#endif

#ifndef CAMERA_H
#define CAMERA_H

#include "vec3.h"
//<position>x y z</position><gaze>x y z</gaze><up>x y z</up><nearplane>left right bottom top</nearplane><neardistance> distance </neardistance><imageresolution>nx ny</imageresolution>

struct Nearplane{
	Nearplane(double left, double right, double top, double bottom): left(left), right(right), top(top), bottom(bottom){};
	double left;
	double right;
	double top;
	double bottom;
};

class Camera{
public:
	Camera(const Point3& position, const Vec3& gaze, const Vec3& up, const Nearplane& nearplane, double neardistance, int nx, int ny): 
		_position(position), _gaze(gaze), _nearplane(nearplane), _neardistance(neardistance), _nx(nx), _ny(ny){
			_w = -unit_vector(gaze);
			_u =  unit_vector(cross(_w, up));
			_v = unit_vector(cross(_w, _u));
		}
private:
	Point3 _position;
	Vec3 _gaze;
	Vec3 _up;
	Nearplane _nearplane;
	double _neardistance;
	int _nx;
	int _ny;
	Vec3 _u;
	Vec3 _v;
	Vec3 _w;

public:
    const Point3& position() const  { return _position; }
    const Vec3& gaze() const  { return _gaze; }
    const Vec3& up() const  { return _up; }
    const Vec3& u() const  { return _u; }
    const Vec3& v() const  { return _v; }
    const Vec3& w() const  { return _w; }
	const Nearplane& nearplane() const {return _nearplane;}
	double neardistance() const{return _neardistance;}
	int nx() const {return _nx;}
	int ny() const {return _ny;}
	void go_u(double amount){_position += amount * _u;}
	void go_v(double amount){_position += amount * _v;}
	void go_w(double amount){_position += amount * _w;}
};

#endif

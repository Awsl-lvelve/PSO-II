#pragma once
#include<random>

class Particle
{
private:
	double _x;
	double _y;
public:
	Particle():_x(0),_y(0){};
	Particle(double x, double y) :_x(x), _y(y) {};

	void set_x(const double& x);
	const double& get_x()const;

	void set_y(const double& y);
	const double& get_y()const;
};
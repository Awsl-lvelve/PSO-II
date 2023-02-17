#include "Particle.h"
#include<time.h>


void Particle::set_x(const double& x)
{
	this->_x = x;
}

const double& Particle::get_x() const
{
	return this->_x;
}

void Particle::set_y(const double& y)
{
	this->_y = y;
}

const double& Particle::get_y() const
{
	return this->_y;
}

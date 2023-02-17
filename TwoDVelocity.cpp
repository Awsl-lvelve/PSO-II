#include"TwoDVelocity.h"

void TwoDVelocity::set_v_x(const double& vx)
{
	this->_v_x = vx;
}

void TwoDVelocity::set_v_y(const double& vy)
{
	this->_v_y = vy;
}

const double& TwoDVelocity::get_v_x() const
{
	return this->_v_x;
}

const double& TwoDVelocity::get_v_y() const
{
	return this->_v_y;
}

double TwoDVelocity::get_combined_v() const
{
	return sqrt(pow(_v_x, 2) + pow(_v_y, 2));
}

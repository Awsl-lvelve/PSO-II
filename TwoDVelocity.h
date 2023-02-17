#pragma once
#include<math.h>

class TwoDVelocity
{
private:
	double _v_x;
	double _v_y;
public:
	TwoDVelocity(const double& vx, const double& vy) :_v_x(vx), _v_y(vy) {};

	TwoDVelocity() :_v_x(1), _v_y(1){};

	void set_v_x(const double&vx);

	void set_v_y(const double& vy);

	const double& get_v_x() const;

	const double& get_v_y()const;

	double get_combined_v()const;
};
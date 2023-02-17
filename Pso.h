#pragma once
#include<iostream>
#include<map>
#include<vector>
#include<fstream>
#include<algorithm>
#include<functional>
#include<random>
#include<cmath>
#include<math.h>

#include"Particle.h"
#include"TwoDVelocity.h"

using namespace std;



class Pso
{
private:
	int _number_of_particles;//������Ŀ
	int _number_of_dimensions;//ά����Ŀ
	int _maximum_of_iteration;//�����������ֵ

	double _x_upper_bound;//�����Ͻ�
	double _x_lower_bound;//�����½�

	double _y_upper_bound;//�����Ͻ�
	double _y_lower_bound;//�����½�


	double _x_velocity_ub;//x�ٶ��Ͻ�
	double _x_velocity_lb;
	double _y_velocity_ub;//�ٶ��½�
	double _y_velocity_lb;

	double _w;//inertia coefficient
	double _c1;//local coefficient
	double _c2;//global coefficient

	vector<Particle>_particles;//����
	vector<TwoDVelocity>_particles_velocity;//�����ٶȼ�¼

	//ȫ������������key���� value.first����ǰȫ��������Ӧ��ֵ��value.second ������ǰȫ�����ɶȵ�����
	map<int, pair<double, Particle>, less<int>>_global_best;


	//�ֲ���������
	//key ���ӱ��
	//value ������ʷ��¼
	//value.key ������Ӧ��ֵ
	//value.value.first ����ֵ
	//value.value.second ��
	map<int, multimap<double, pair<Particle, int>, less<double>>, less<int>>_individual_best;


public:

	friend ostream& operator<<(ostream& cout, Particle p);

	Pso(int NOP = 50,
		int NOD = 1,
		int MOI = 1000,
		double XUB = 5.0, double XLB = -5.0,
		double YUB = 5.0, double YLB = -5.0,
		double VXUB = 0.5, double VXLB = -0.5,
		double VYUB = 0.5, double VYLB = -0.5,
		double W = 0.4, 
		double C1 = 0.5, 
		double C2 = 0.2) :

		_number_of_particles(NOP),
		_number_of_dimensions(NOD),
		_maximum_of_iteration(MOI),
		
		_x_upper_bound(XUB),_x_lower_bound(XLB),
		_y_upper_bound(YUB),_y_lower_bound(YLB),
		_x_velocity_ub(VXUB),_x_velocity_lb(VXLB),
		_y_velocity_ub(VYUB),_y_velocity_lb(VYLB),
		
		_w(W),
		_c1(C1),
		_c2(C2)
	{
		this->_particles.resize(NOP, Particle());
		this->_particles_velocity.resize(NOP, TwoDVelocity());
	}

	void set_nop(const int& nop);
	void set_nod(const int& nod);
	void set_moi(const int& moi);


	void set_w(const double& w);
	void set_c1(const double& c1);
	void set_c2(const double& c2);

	//��ʼ������


	void calc_v();

	const double& calc_fit_value(Particle p);

	void init();

	void particle_init();

	void individual_best_init();

	void global_best_init();

	void evolve();

	void optimize();

	void update_individual_best();

	void update_global_best();

	void print_particles();

	void print_individual_best();

	void print_global_best();

	void print_individual_results();

	void print_global_result();

	void save_global_best();

	void save_individual_best();
};
#include"Pso.h"

void Pso::set_nop(const int& nop)
{
	this->_number_of_particles = nop;
}

void Pso::set_nod(const int& nod)
{
	this->_number_of_dimensions = nod;
}

void Pso::set_moi(const int& moi)
{
	this->_maximum_of_iteration = moi;
}





void Pso::set_w(const double& w)
{
	this->_w = w;
}

void Pso::set_c1(const double& c1)
{
	this->_c1 = c1;
}

void Pso::set_c2(const double& c2)
{
	this->_c2 = c2;
}



void Pso::calc_v()
{
	map<int, multimap<double, pair<Particle, int>, less<double>>, less<int>>::iterator iiter = this->_individual_best.begin();
	
	auto giter = this->_global_best.rbegin();

	vector<Particle>::iterator piter = this->_particles.begin();

	for (vector<TwoDVelocity>::iterator iter = this->_particles_velocity.begin();

		iter != this->_particles_velocity.end() &&
		piter != this->_particles.end() &&
		iiter != this->_individual_best.end();

		iter++,
		piter++,
		iiter++) {

		mt19937 mersenne(static_cast<unsigned int>(time(nullptr)));
		uniform_real_distribution<double>rnd_global(0, 1);

		double c_arr[4] = { 0.0, 0.0,0.0,0.0 };

		for (int i = 0; i < 4; i++)
		{
			c_arr[i] = rnd_global(mersenne);
		}

		double v_x = iter->get_v_x() * this->_w +
			this->_c1 * c_arr[0] * (iiter->second.begin()->first - piter->get_x()) +
			this->_c2 * c_arr[1] * (giter->second.first - piter->get_x());

		double v_y = iter->get_v_x() * this->_w +
			this->_c1 * c_arr[2] * (iiter->second.begin()->first - piter->get_y()) +
			this->_c2 * c_arr[3] * (giter->second.first - piter->get_y());

		iter->set_v_x(v_x);
		iter->set_v_y(v_y);

		if (iter->get_v_x()>this->_x_velocity_ub)
		{
			iter->set_v_x(this->_x_velocity_ub);
		}
		else if(iter->get_v_x()<this->_x_velocity_lb)
		{
			iter->set_v_x(this->_x_velocity_lb);
		}

		if (iter->get_v_y()>this->_y_velocity_ub)
		{
			iter->set_v_y(this->_y_velocity_ub);
		}
		else if (iter->get_v_y()<_y_velocity_lb)
		{
			iter->set_v_y(this->_y_velocity_lb);
		}

	}
}

const double& Pso::calc_fit_value(Particle p)
{
	double val = 20 + pow(p.get_x(), 2) + pow(p.get_y(), 2)
		- 10 * cos(2 * 3.141592653579 * p.get_x())
		- 10 * cos(2 * 3.141592653579 * p.get_y());

	return val;
}

void Pso::init()
{
	this->particle_init();
	this->individual_best_init();
	this->global_best_init();


}

void Pso::particle_init()
{
	cout << "Particles are initing..." << endl;

	mt19937 mersenne(static_cast<unsigned int>(time(nullptr)));
	for (vector<Particle>::iterator iter = this->_particles.begin();
		iter != this->_particles.end();
		iter++) {
		
	
		uniform_real_distribution<double>rnd_x(this->_x_lower_bound,this->_x_upper_bound);//产生均匀分布
		double rnd_num_x = rnd_x(mersenne);

		uniform_real_distribution<double>rnd_y(this->_y_lower_bound, this->_y_upper_bound);
		double rnd_num_y = rnd_y(mersenne);
			
		iter->set_x(rnd_num_x);
		iter->set_y(rnd_num_y);
	}

}

void Pso::individual_best_init()
{
	this->_individual_best;
	int i = 0;
	for (vector<Particle>::iterator iter = this->_particles.begin();
		iter != this->_particles.end(); iter++,i++) {
		multimap<double, pair<Particle, int>, less<double>>m;

		double fit_val = calc_fit_value(*iter);

		m.insert(make_pair(fit_val, pair<Particle, int>(*iter, 1)));

		this->_individual_best.insert(make_pair(i, m));
	}
}

void Pso::global_best_init()
{
	mt19937 mersenne(static_cast<unsigned int>(time(nullptr)));
	uniform_real_distribution<double>rnd_global(0, 1);
	this->_global_best;

	double x = this->_x_lower_bound + (this->_x_upper_bound - this->_x_lower_bound) * rnd_global(mersenne);
	double y = this->_y_lower_bound + (this->_y_upper_bound - this->_y_lower_bound) * rnd_global(mersenne);

	Particle particle(x, y);

	
	double test_g=this->calc_fit_value(particle);

	pair<double, Particle>p(test_g,particle);


	this->_global_best.insert(make_pair(1, p));
}

ostream& operator<<(ostream& cout, Particle p)
{
	cout << "(" << p.get_x() << "," << p.get_y() << ")";
	return cout;
}

void Pso::evolve()
{
	this->calc_v();
	vector<TwoDVelocity>::iterator viter = this->_particles_velocity.begin();

	for (vector<Particle>::iterator piter = this->_particles.begin(); piter != this->_particles.end() && viter != this->_particles_velocity.end(); piter++, viter++) {
		
		piter->set_x(piter->get_x() + viter->get_v_x());
		piter->set_y(piter->get_y() + viter->get_v_y());

		if (piter->get_x()>this->_x_upper_bound)
		{
			piter->set_x(this->_x_upper_bound);
		}
		else if(piter->get_x()<this->_x_lower_bound)
		{
			piter->set_x(this->_x_lower_bound);
		}

		if (piter->get_y()>this->_y_upper_bound)
		{
			piter->set_y(this->_y_upper_bound);
		}
		else if(piter->get_y()<this->_y_lower_bound)
		{
			piter->set_y(this->_y_lower_bound);
		}
		
	}
}

void Pso::optimize()
{
	int i = 1;
	while (i<=this->_maximum_of_iteration)
	{
		evolve();
		update_individual_best();
		update_global_best();
		i++;
	}
}

void Pso::update_individual_best()
{
	map<int, multimap<double, pair<Particle, int>, less<double>>, less<int>>::iterator miter = this->_individual_best.begin();

	for (vector<Particle>::iterator iter = this->_particles.begin(); 
		iter != this->_particles.end()&&miter!=this->_individual_best.end();
		miter++,iter++) {
		double fit_val = calc_fit_value(*iter);

		pair<Particle, int>p = miter->second.begin()->second;
		p.second = miter->second.size() + 1;//轮数等于粒子记录的条数++1 每轮都会更新

		double current_best_fit_val = miter->second.begin()->first;//当前最优值

		if (fit_val < miter->second.begin()->first)//如果新的粒子的适应度值更小
		{
			p.first = *iter;//更新最优值

			miter->second.insert(make_pair(fit_val, p));//插入map
		}
		else
		{
			miter->second.insert(make_pair(current_best_fit_val, p));//轮数标记加1，但是最优值和适应度值不发生变化
		}

	}

}

void Pso::update_global_best()
{
	pair<double, Particle>current_global_best = this->_global_best.rbegin()->second;

	map<int, multimap<double, pair<Particle, int>, less<double>>, less<int>>::iterator iter;

	for (iter = this->_individual_best.begin(); iter != this->_individual_best.end(); iter++) {
		if (current_global_best.first > iter->second.begin()->first)
		{
			current_global_best.first = iter->second.begin()->first;//更新适应度
			current_global_best.second = iter->second.begin()->second.first;//更新值
		}

	}
	this->_global_best.insert(make_pair(this->_global_best.size() + 1, current_global_best));
}

void Pso::print_particles()
{
	for (auto& particle : this->_particles) {
		cout << particle << " ";
	}
	cout << endl;
}

void Pso::print_individual_best()
{
	cout << "粒子编号 " << "粒子适应度 " << "粒子值 " << "轮" << endl;
	for (map<int, multimap<double, pair<Particle, int>, less<double>>, less<int>>::iterator iter = this->_individual_best.begin();
		iter != this->_individual_best.end();
		iter++) {
		for (multimap<double, pair<Particle, int>, less<double>>::iterator iter2 = iter->second.begin();
			iter2 != iter->second.end();
			iter2++) {
			cout << iter->first << " ";
			cout << iter2->first << " ";
			cout << iter2->second.first << " ";
			cout << iter2->second.second << endl;
		}
	}
}

void Pso::print_global_best(){

	cout << "全局最优" << endl;
	cout << "轮 " << "粒子适应度 " << "粒子值" << endl;
	for (auto& c : this->_global_best) {
		cout << c.first << " ";
		cout << c.second.first << " ";
		cout << c.second.second << endl;
	}
}

void Pso::print_individual_results()
{
}

void Pso::print_global_result()
{
}

void Pso::save_global_best()
{
	ofstream ofs("global_best_recs.txt", ios::out | ios::trunc);

	for (auto& c : this->_global_best) {
		ofs << c.first << " ";
		ofs << c.second.first << " ";
		ofs << c.second.second << endl;
	}

	ofs.close();
}

void Pso::save_individual_best(){

	ofstream ofs("individual_best.txt", ios::out | ios::trunc);
	for (map<int, multimap<double, pair<Particle, int>, less<double>>, less<int>>::iterator particle_it = this->_individual_best.begin(); particle_it != this->_individual_best.end(); particle_it++) {
		for (auto& deatil : particle_it->second) {
			ofs << particle_it->first << ' ';
			ofs << deatil.first << ' ';
			ofs << deatil.second.first << ' ';
			ofs << deatil.second.second << endl;
		}
	}
	ofs.close();


}

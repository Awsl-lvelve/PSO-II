#include"Pso.h"

int main() {
	Pso pso;
	pso.init();
	pso.print_particles();
	pso.print_individual_best();
	pso.print_global_best();

	pso.optimize();

	pso.save_global_best();
	pso.save_individual_best();

	pso.print_global_result();
}
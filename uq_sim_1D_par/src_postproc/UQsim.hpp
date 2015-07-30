#ifndef UQSIM_HPP_
#define UQSIM_HPP_

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <cstdlib>

class UQSimulation
{
protected:
	int rank;
	int nprocs;

	double rho_f_p1;
	double rho_f_p2;
	double nu_f_p1;
	double nu_f_p2;
	double rho_s_p1;
	double rho_s_p2;
	std::string nastin_dat;
	std::string solidz_dat;
	std::string create_data_rank;
	std::string run_exec;
	std::string output_data;
	std::string gather_data_exec_mc;
	std::string postproc_stat_exec_mc;
	std::string postproc_file_all_mc;
	std::string postproc_stat_mc;
	std::string gather_data_exec_sc;
	std::string get_output_sc;
	std::string postproc_stat_exec_sc;
	std::string output_file_sc;
	std::string coeff_sc;
	std::string postproc_stat_sc;
	std::string insert_nastin_exec;
	std::string insert_solidz_exec;
	std::string gather_alya_output;

public:
	virtual int local_global_mapping(const int local_index, const int& rank) const = 0;
	virtual int data_decomp() const = 0;

	virtual std::vector<double> pre_processing() const = 0;
	virtual std::vector<double> pre_processing(const double& param1, const double& param2) const = 0;
	
	virtual void simulation(std::vector<double>& pre_proc_result) const = 0;

	virtual void post_processing() const = 0;

	virtual ~UQSimulation() {}
};

#endif /* UQSIM_HPP_ */
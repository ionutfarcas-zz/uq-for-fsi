#ifndef MCSIMNORMAL_HPP_
#define MCSIMNORMAL_HPP_

#include "UQsim.hpp"
#include "normalRV.hpp"
#include "helper.hpp"

class MCSimulation_normal : public UQSimulation
{
private:
	int nsamples;
	int local_samples;
	double mean;
	double std_dev;
	
	std::string create_data_each_rank;
	int run_create_data_rank_ok;

	NormalRandomVariable nrv;

public:
	MCSimulation_normal() 
	{
		nsamples = 0;
		local_samples = 0;
		mean = 0.0;
		std_dev = 0.0;
	}

	MCSimulation_normal(
		std::string& _nastin_dat, 
		std::string& _solidz_dat,
		std::string& _create_data_rank, 
		std::string& _run_exec,
		std::string& _output_data, 
		std::string& _gather_data_exec_mc, 
		std::string& _postproc_stat_exec_mc,  
		std::string& _postproc_file_all_mc, 
		std::string& _postproc_stat_mc, 
		std::string& _insert_nastin_exec, 
		std::string& _insert_solidz_exec,
		std::string& _gather_alya_output, 
		const unsigned int& _nsamples,
		const unsigned int& _rank,
		const unsigned int& _nprocs, 
		const double& _rho_f_p1,  
		const double& _rho_f_p2, 
		const double& _nu_f_p1, 
		const double& _nu_f_p2, 
		const double& _rho_s_p1, 
		const double& _rho_s_p2) 
	{
		mean = 0.0;
		std_dev = 1.0;
		
		nastin_dat = _nastin_dat;
		solidz_dat = _solidz_dat;
		create_data_rank = _create_data_rank;
		run_exec = _run_exec;
		output_data = _output_data;
		gather_data_exec_mc = _gather_data_exec_mc;
		postproc_stat_exec_mc = _postproc_stat_exec_mc;
		postproc_file_all_mc = _postproc_file_all_mc;
		postproc_stat_mc = _postproc_stat_mc;
		insert_nastin_exec = _insert_nastin_exec;
		insert_solidz_exec = _insert_solidz_exec;
		gather_alya_output = _gather_alya_output;

		assert(_nsamples >= _nprocs);
		nsamples = _nsamples;
		rank = _rank;
		nprocs = _nprocs;

		local_samples = data_decomp();

		rho_f_p1 = _rho_f_p1;
		rho_f_p2 = _rho_f_p2;
		nu_f_p1 = _nu_f_p1;
		nu_f_p2 = _nu_f_p2;
		rho_s_p1 = _rho_s_p1;
		rho_s_p2 = _rho_s_p2;

		create_data_each_rank = run_create_data_rank(create_data_rank, rank);
		run_create_data_rank_ok = system(create_data_each_rank.c_str());
		assert(run_create_data_rank_ok >=0 );
	}

	virtual int local_global_mapping(const int local_index, const int& rank) const
	{
		int global_id;
		int local_samples;

		local_samples = nsamples/nprocs;
		global_id = local_index + rank*local_samples;

		return global_id;
	}

	virtual int data_decomp() const
	{
		int local_points = 0;
		int rem = 0;

		local_points = nsamples/nprocs;
		rem = nsamples % nprocs;

		if (rank == nprocs - 1)   
		{
			local_points += rem;
		}

		return local_points;
	}

	virtual std::vector<double> pre_processing() const
	{
		std::vector<double> samples;

		samples = nrv.get_samples(mean, std_dev, local_samples);

		return samples;
	}

	virtual std::vector<double> pre_processing(const double& param1, const double& param2) const
	{
		std::vector<double> samples;

		samples = nrv.get_samples(param1, param2, local_samples);

		return samples;
	}

	virtual void simulation(std::vector<double>& pre_proc_result) const  
	{
		std::string modify_nastin_data;
		int modify_nastin_data_ok = 0;

		double rand_par = 0.0;

		for(int i = 0 ; i < local_samples ; ++i)
		{
			rand_par = pre_proc_result[i];

			modify_nastin_data = run_insert_nastin_1d(insert_nastin_exec, nastin_dat, rand_par, rank);
			modify_nastin_data_ok = system(modify_nastin_data.c_str());
			assert(modify_nastin_data_ok >= 0);
		}
	}

	virtual void post_processing() const
	{
		
	}

	~MCSimulation_normal() {}
};

#endif /* MCSIMNORMAL_HPP_ */
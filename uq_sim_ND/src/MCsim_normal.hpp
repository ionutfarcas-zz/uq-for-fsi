#ifndef MCSIMNORMAL_HPP_
#define MCSIMNORMAL_HPP_

#include "UQsim.hpp"
#include "normalRV.hpp"

class MCSimulation_normal : public UQSimulation
{
private:
	int nsamples;
	double mean;
	double std_dev;

	NormalRandomVariable nrv;

public:
	MCSimulation_normal() 
	{
		nsamples = 0;
		mean = 0.0;
		std_dev = 0.0;
	}

	MCSimulation_normal(
		std::string& _nastin_dat, 
		std::string& _solidz_dat, 
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
		run_exec = _run_exec;
		output_data = _output_data;
		gather_data_exec_mc = _gather_data_exec_mc;
		postproc_stat_exec_mc = _postproc_stat_exec_mc;
		postproc_file_all_mc = _postproc_file_all_mc;
		postproc_stat_mc = _postproc_stat_mc;
		insert_nastin_exec = _insert_nastin_exec;
		insert_solidz_exec = _insert_solidz_exec;
		gather_alya_output = _gather_alya_output;

		nsamples = _nsamples;

		rho_f_p1 = _rho_f_p1;
		rho_f_p2 = _rho_f_p2;
		nu_f_p1 = _nu_f_p1;
		nu_f_p2 = _nu_f_p2;
		rho_s_p1 = _rho_s_p1;
		rho_s_p2 = _rho_s_p2;
	}
	
	virtual double compute_volume() const
	{
		double volume = 1.0;

		return volume;
	}

	virtual vec2d_float_t pre_processing() const
	{
		vec2d_float_t vec;

		return vec;
	}

	virtual std::vector<double> pre_processing(const double& param1, const double& param2) const
	{
		std::vector<double> samples;

		samples = nrv.get_samples(param1, param2, nsamples);

		return samples;
	}

	virtual void simulation(const vec2d_float_t& pre_proc_result) const
	{
		/*TO DO - Nothing here */
	}

	virtual void simulation(std::vector<double>& pre_proc_result_dim1, std::vector<double>& pre_proc_result_dim2) const
	{
		std::string modify_nastin_data;
		std::string modify_solidz_data;
		std::string get_alya_output;
		std::string get_all_data;

		int modify_nastin_data_ok = 0;
		int modify_solidz_data_ok = 0;
		int run_ok = 0;
		int gather_alya_output_ok = 0;
		int get_all_data_ok = 0;

		for(int i = 0 ; i < nsamples ; ++i)
		{
			modify_nastin_data = run_insert_nastin_1d(insert_nastin_exec, nastin_dat, pre_proc_result_dim1[i]);
			modify_nastin_data_ok = system(modify_nastin_data.c_str());
			assert(modify_nastin_data_ok >= 0);

			modify_solidz_data = run_insert_solidz_1d(insert_solidz_exec, solidz_dat, pre_proc_result_dim2[i]);
			modify_solidz_data_ok= system(modify_solidz_data.c_str());
			assert(modify_solidz_data_ok >= 0);

			run_ok = system(run_exec.c_str());
			assert(run_ok >= 0);

			get_alya_output = run_gather_alya_output(gather_alya_output, i+1);
			gather_alya_output_ok = system(get_alya_output.c_str());
			assert(gather_alya_output_ok >=0);

			get_all_data = run_gather_data(gather_data_exec_mc, output_data, postproc_file_all_mc, i+1);
			get_all_data_ok = system(get_all_data.c_str());
			assert(get_all_data_ok >= 0);
		}
	}

	virtual void post_processing() const
	{
		std::string get_postproc_stat;
		int get_postproc_stat_ok = 0;

		get_postproc_stat = run_postproc_stat(postproc_stat_exec_mc, postproc_file_all_mc, postproc_stat_mc, nsamples);
		get_postproc_stat_ok = system(get_postproc_stat.c_str());
		assert(get_postproc_stat_ok >= 0);
	}

	~MCSimulation_normal() {}
};

#endif /* MCSIMNORMAL_HPP_ */
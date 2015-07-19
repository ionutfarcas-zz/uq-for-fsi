#ifndef MCSIMUNIFORM_HPP_
#define MCSIMUNIFORM_HPP_

#include "UQsim.hpp"
#include "uniformRV.hpp"
#include "helper.hpp"

class MCSimulation_uniform : public UQSimulation
{
private:
	int nsamples;
	int local_samples;
	double left_param;
	double right_param;

	std::string create_data_each_rank;
	int run_create_data_rank_ok;

	UniformRandomVariable urv;

public:
	MCSimulation_uniform() 
	{
		nsamples = 0;
		local_samples = 0;
		left_param = 0.0;
		right_param = 0.0;
	}

	MCSimulation_uniform(
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
		left_param = 0.0;
		right_param = 1.0;
		
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

	virtual double compute_volume() const
	{
		return 1.0;
	}

	virtual vec2d_float_t pre_processing() const
	{
		vec2d_float_t samples;

		return samples;
	}

	virtual std::vector<double> pre_processing(const double& param1, const double& param2) const
	{
		std::vector<double> samples;

		samples = urv.get_samples(param1, param2, local_samples);

		return samples;
	}

	virtual void simulation(const vec2d_float_t& pre_proc_result) const 
	{
		/*TO DO : nothing here */
	}

	virtual void simulation(std::vector<double>& pre_proc_result_dim1, std::vector<double>& pre_proc_result_dim2) const
	{
		int global_id = 0;

		std::string modify_nastin_data;
		std::string modify_solidz_data;
		std::string alya_nastin_solidz;
		std::string get_alya_output;
		std::string get_all_data;

		int modify_nastin_data_ok = 0;
		int modify_solidz_data_ok = 0;
		int run_alya_ok = 0;
		int gather_alya_output_ok = 0;
		int get_all_data_ok = 0;

		double rand_par_dim1 = 0.0;
		double rand_par_dim2 = 0.0;

		for(int i = 0 ; i < local_samples ; ++i)
		{
			global_id = this->local_global_mapping(i, rank);
			rand_par_dim1 = pre_proc_result_dim1[i];
			rand_par_dim2 = pre_proc_result_dim2[i];

			modify_nastin_data = run_insert_nastin_1d(insert_nastin_exec, nastin_dat, rand_par_dim1, rank);
			modify_nastin_data_ok = system(modify_nastin_data.c_str());
			assert(modify_nastin_data_ok >= 0);
			modify_solidz_data = run_insert_solidz_1d(insert_solidz_exec, solidz_dat, rand_par_dim2, rank);
			modify_solidz_data_ok= system(modify_solidz_data.c_str());
			assert(modify_solidz_data_ok >= 0);

			alya_nastin_solidz = run_alya(run_exec, rank);
			run_alya_ok = system(alya_nastin_solidz.c_str());
			assert(run_alya_ok >= 0);

			get_alya_output = run_gather_alya_output(gather_alya_output, global_id+1);
			gather_alya_output_ok = system(get_alya_output.c_str());
			assert(gather_alya_output_ok >=0);

			get_all_data = run_gather_data(gather_data_exec_mc, output_data, postproc_file_all_mc, global_id+1, rank);
			get_all_data_ok = system(get_all_data.c_str());
			assert(get_all_data_ok >= 0);
		}
	}

	virtual void post_processing() const
	{
		std::vector<double> sums;

		double temp_sum_disp_x = 0.0;
		double temp_sum_forces_x = 0.0;
		double temp_sum_forces_y = 0.0;
		double temp_sum2_disp_x = 0.0;
		double temp_sum2_forces_x = 0.0;
		double temp_sum2_forces_y = 0.0;

		double sum_disp_x = 0.0;
		double sum_forces_x = 0.0;
		double sum_forces_y = 0.0;
		double sum2_disp_x = 0.0;
		double sum2_forces_x = 0.0;
		double sum2_forces_y = 0.0;
		int temp_no_valid_lines = 0;
		int no_valid_lines = 0;

		double mean_disp_x = 0.0;
		double mean_forces_x = 0.0;
		double mean_forces_y = 0.0;
		double var_disp_x = 0.0;
		double var_forces_x = 0.0;
		double var_forces_y = 0.0;

		std::string get_postproc_stat;

		get_postproc_stat = run_postproc_stat(postproc_stat_exec_mc, postproc_file_all_mc, rank);
		
		sums = get_stat_mc(get_postproc_stat, temp_no_valid_lines);
		temp_sum_disp_x = sums[0];
		temp_sum_forces_x = sums[1];
		temp_sum_forces_y = sums[2];
		temp_sum2_disp_x = sums[3];
		temp_sum2_forces_x = sums[4];
		temp_sum2_forces_y = sums[5];

		MPI_Allreduce(&temp_no_valid_lines, &no_valid_lines, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		MPI_Allreduce(&temp_sum_disp_x, &sum_disp_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&temp_sum_forces_x, &sum_forces_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&temp_sum_forces_y, &sum_forces_y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&temp_sum2_disp_x, &sum2_disp_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&temp_sum2_forces_x, &sum2_forces_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&temp_sum2_forces_y, &sum2_forces_y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		mean_disp_x = sum_disp_x/no_valid_lines;
		mean_forces_x = sum_forces_x/no_valid_lines;
		mean_forces_y = sum_forces_y/no_valid_lines;

		var_disp_x = (sum2_disp_x - (sum_disp_x*sum_disp_x)/no_valid_lines)/(no_valid_lines - 1.0);
		var_forces_x = (sum2_forces_x - (sum_forces_x*sum_forces_x)/no_valid_lines)/(no_valid_lines - 1.0);
		var_forces_y = (sum2_forces_y - (sum_forces_y*sum_forces_y)/no_valid_lines)/(no_valid_lines - 1.0);

		if(rank == 0)
		{
			int save_stats_ok = 0;
			save_stats_ok = write_stat_to_file(postproc_stat_mc, mean_disp_x, mean_forces_x, mean_forces_y, var_disp_x, var_forces_x, var_forces_y);
			assert(save_stats_ok == 0);
		}
	}
	
	~MCSimulation_uniform() {}
};

#endif /* MCSIMUNIFORM_HPP_ */
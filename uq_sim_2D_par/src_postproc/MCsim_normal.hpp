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
		/*TO DO : nothing here */
	}

	virtual void simulation(std::vector<double>& pre_proc_result_dim1, std::vector<double>& pre_proc_result_dim2) const
	{
		int global_id = 0;

		std::string get_alya_output;
		std::string get_all_data;

		int gather_alya_output_ok = 0;
		int get_all_data_ok = 0;

		for(int i = 0 ; i < local_samples ; ++i)
		{
			global_id = this->local_global_mapping(i, rank);	

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

		double mean_disp_x = 0.0;
		double mean_forces_x = 0.0;
		double mean_forces_y = 0.0;
		double var_disp_x = 0.0;
		double var_forces_x = 0.0;
		double var_forces_y = 0.0;

		int local_number_of_datapoints = 0;
		int local_time_steps = 0;

		std::string get_postproc_stat;

		get_postproc_stat = run_postproc_stat(postproc_stat_exec_mc, postproc_file_all_mc, local_samples, rank);
		
		sums = get_stat_mc(get_postproc_stat, local_number_of_datapoints);
		local_time_steps = local_number_of_datapoints/local_samples;

		for(int i = 0 ; i < local_time_steps ; ++i)
		{
			temp_sum_disp_x = sums[i*6];
			temp_sum_forces_x = sums[i*6 + 1];
			temp_sum_forces_y = sums[i*6 + 2];
			temp_sum2_disp_x = sums[i*6 + 3];
			temp_sum2_forces_x = sums[i*6 + 4];
			temp_sum2_forces_y = sums[i*6 + 5];

			MPI_Allreduce(&temp_sum_disp_x, &sum_disp_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&temp_sum_forces_x, &sum_forces_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&temp_sum_forces_y, &sum_forces_y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&temp_sum2_disp_x, &sum2_disp_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&temp_sum2_forces_x, &sum2_forces_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&temp_sum2_forces_y, &sum2_forces_y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			mean_disp_x = sum_disp_x/nsamples;
			mean_forces_x = sum_forces_x/nsamples;
			mean_forces_y = sum_forces_y/nsamples;

			var_disp_x = (sum2_disp_x - (sum_disp_x*sum_disp_x)/nsamples)/(nsamples - 1.0);
			var_forces_x = (sum2_forces_x - (sum_forces_x*sum_forces_x)/nsamples)/(nsamples - 1.0);
			var_forces_y = (sum2_forces_y - (sum_forces_y*sum_forces_y)/nsamples)/(nsamples - 1.0);

			if(rank == 0)
			{
				int save_stats_ok = 0;
				int time_step = i + 1;

				save_stats_ok = write_stat_to_file(postproc_stat_mc, mean_disp_x, mean_forces_x, mean_forces_y, var_disp_x, var_forces_x, var_forces_y, time_step);
				assert(save_stats_ok == 0);
			}
		}
	}

	~MCSimulation_normal() {}
};

#endif /* MCSIMNORMAL_HPP_ */
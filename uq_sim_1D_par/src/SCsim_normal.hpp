#ifndef SCSIMNORMAL_HPP_
#define SCSIMNORMAL_HPP_

#include "UQsim.hpp"
#include "gauss_hermite_quad.hpp"
#include "helper.hpp"

class SCSimulation_normal : public UQSimulation
{
private:
	int ncoeff;
	int quad_degree;
	int local_quad_points;

	std::string create_data_each_rank;
	int run_create_data_rank_ok;
	
	GaussHermiteQuadrature ghq;

public:
	SCSimulation_normal() 
	{
		ncoeff = 0;
		quad_degree = 0;
		local_quad_points = 0;

		create_data_each_rank = " ";
		run_create_data_rank_ok = 0;
	}

	SCSimulation_normal(
		std::string& _nastin_dat, 
		std::string& _solidz_dat,
		std::string& _create_data_rank,  
		std::string& _run_exec, 
		std::string& _output_data, 
		std::string& _gather_data_exec_sc, 
		std::string& _get_output_sc,
		std::string& _postproc_stat_exec_sc,
		std::string& _output_file_sc, 
		std::string& _coeff_sc,
		std::string& _postproc_stat_sc,
		std::string& _insert_nastin_exec, 
		std::string& _insert_solidz_exec,
		std::string& _gather_alya_output,  
		const unsigned int& _ncoeff, 
		const unsigned int& _quad_degree,
		unsigned int _rank,
		const unsigned int& _nprocs,
		const double& _rho_f_p1, 
		const double& _rho_f_p2, 
		const double& _nu_f_p1, 
		const double& _nu_f_p2, 
		const double& _rho_s_p1, 
		const double& _rho_s_p2)
	{
		ncoeff = _ncoeff;
		quad_degree = _quad_degree;
		
		nastin_dat = _nastin_dat;
		solidz_dat = _solidz_dat;
		create_data_rank = _create_data_rank;
		run_exec = _run_exec;
		output_data = _output_data;
		gather_data_exec_sc = _gather_data_exec_sc;
		get_output_sc = _get_output_sc;
		postproc_stat_exec_sc = _postproc_stat_exec_sc;
		output_file_sc = _output_file_sc;
		coeff_sc = _coeff_sc;
		postproc_stat_sc = _postproc_stat_sc;
		insert_nastin_exec = _insert_nastin_exec;
		insert_solidz_exec = _insert_solidz_exec;
		gather_alya_output = _gather_alya_output;

		assert(_quad_degree >= _nprocs);
		rank = _rank;
		nprocs = _nprocs;

		create_data_each_rank = run_create_data_rank(create_data_rank, _rank);
		run_create_data_rank_ok = system(create_data_each_rank.c_str());
		assert(run_create_data_rank_ok >=0 );

		local_quad_points = data_decomp();

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
		int local_quad_degree;

		local_quad_degree = quad_degree/nprocs;
		global_id = local_index + rank*local_quad_degree;

		return global_id;
	}

	virtual int data_decomp() const
	{
		int local_quad_degree = 0;
		int rem = 0;

		local_quad_degree = quad_degree/nprocs;
		rem = quad_degree % nprocs;

		if (rank == nprocs - 1)   
		{
			local_quad_degree += rem;
		}

		return local_quad_degree;
	} 

	virtual std::vector<double> pre_processing() const
	{
		std::vector<double> nodes, weights, nodes_weights;

		ghq.quad_nodes_weights(quad_degree, weights, nodes);

		for(int i = 0 ; i < quad_degree ; ++i)
		{
			nodes_weights.push_back(nodes[i]);
		}

		for(int i = 0 ; i < quad_degree ; ++i)
		{
			nodes_weights.push_back(weights[i]);
		}

		return nodes_weights;
	}

	virtual std::vector<double> pre_processing(const double& param1, const double& param2) const
	{
		// do nothing; useful only for MCS
		std::vector<double> dummy;

		return dummy;
	}
	
	virtual void simulation(std::vector<double>& pre_proc_result) const
	{
		int global_id = 0;

		std::string modify_nastin_data;
		std::string alya_nastin_solidz;
		std::string get_data;
		std::string get_alya_output;
		std::string get_output;

		int run_alya_ok = 0;
		int modify_nastin_data_ok = 0;
		int get_data_ok = 0;
		int gather_alya_output_ok = 0;

		int no_valid_lines_local = 0;
		int no_valid_lines = 0;

		double temp = 0.0;
		double temp_disp_x = 0.0;
		double temp_force0 = 0.0;
		double temp_force1 = 0.0;
		double disp_x_local = 0.0;
		double force0_local = 0.0;
		double force1_local = 0.0;
		double disp_x = 0.0;
		double force0 = 0.0;
		double force1 = 0.0;

		std::vector<double> disp_x_all;
		std::vector<double> force0_all;
		std::vector<double> force1_all;
		vec2d_double get_output_values;

		for (int i = 0; i < local_quad_points; ++i)
		{
			global_id = this->local_global_mapping(i, rank);

			temp = sqrt(2.0)*pre_proc_result[global_id]*nu_f_p2 + nu_f_p1;
			assert(temp >= 0);
			modify_nastin_data = run_insert_nastin_1d(insert_nastin_exec, nastin_dat, temp, rank);
			modify_nastin_data_ok = system(modify_nastin_data.c_str());
			assert(modify_nastin_data_ok >= 0);

			alya_nastin_solidz = run_alya(run_exec, rank);
			run_alya_ok = system(alya_nastin_solidz.c_str());
			assert(run_alya_ok >= 0);

			get_alya_output = run_gather_alya_output(gather_alya_output, global_id+1);
			gather_alya_output_ok = system(get_alya_output.c_str());
			assert(gather_alya_output_ok >=0);

			get_data = run_gather_data(gather_data_exec_sc, output_data, output_file_sc, global_id+1, rank);
			get_data_ok = system(get_data.c_str());
			assert(get_data_ok >= 0);

			get_output = run_get_output(get_output_sc, output_file_sc, rank);
			get_output_values = get_output_data(get_output, no_valid_lines_local);

			disp_x_local = get_output_values[i][0];
			force0_local = get_output_values[i][1];
			force1_local = get_output_values[i][2];

			disp_x_all.push_back(disp_x_local);
			force0_all.push_back(force0_local);
			force1_all.push_back(force1_local);
		}

		MPI_Allreduce(&no_valid_lines_local, &no_valid_lines, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		/* Insert fault-tolerant mechanism here by reducing the quadrature degre */

		for(int j = 0 ; j < ncoeff ; ++j)
		{
			temp_disp_x = 0.0;
			temp_force0 = 0.0;
			temp_force1 = 0.0;

			for (int i = 0 ; i < local_quad_points ; ++i)
			{
				temp_disp_x += disp_x_all[i]*ghq.orthogonal_poly(j, sqrt(2.0)*pre_proc_result[global_id]) * pre_proc_result[quad_degree + global_id]/sqrt(M_PI);
				temp_force0 += force0_all[i]*ghq.orthogonal_poly(j, sqrt(2.0)*pre_proc_result[global_id]) * pre_proc_result[quad_degree + global_id]/sqrt(M_PI);
				temp_force1 += force1_all[i]*ghq.orthogonal_poly(j, sqrt(2.0)*pre_proc_result[global_id]) * pre_proc_result[quad_degree + global_id]/sqrt(M_PI);
			}

			temp_disp_x = temp_disp_x/ghq.norm_factor(j);
			temp_force0 = temp_force0/ghq.norm_factor(j);
			temp_force1 = temp_force1/ghq.norm_factor(j);

			MPI_Allreduce(&temp_disp_x, &disp_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&temp_force0, &force0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&temp_force1, &force1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			if(rank == 0)
			{
				int save_coeff_ok = save_coeff(coeff_sc, disp_x, force0, force1);
				assert(save_coeff_ok == 1);	
			}	
		}
	}

	virtual void post_processing() const
	{
		int get_coeff_ok = 0;
		std::vector<double> disp_x;
		std::vector<double> force0;
		std::vector<double> force1;

		double mean_disp_x = 0.0;
		double mean_forces_x = 0.0;
		double mean_forces_y = 0.0;
		double var_disp_x = 0.0;
		double var_forces_x = 0.0;
		double var_forces_y = 0.0;

		get_coeff_ok = get_coeff_sc(coeff_sc, disp_x, force0, force1);
		assert(get_coeff_ok == 1);

		mean_disp_x = disp_x[0];
		mean_forces_x = force0[0];
		mean_forces_y = force1[0];

		for(int i = 1 ; i < ncoeff ; ++i)
		{
			var_disp_x += disp_x[i]*disp_x[i];
			var_forces_x += force0[i]*force0[i];
			var_forces_y += force1[i]*force1[i];
		}

		if(rank == 0)
		{
			int save_stats_ok = 0;
			save_stats_ok = write_stat_to_file(postproc_stat_sc, mean_disp_x, mean_forces_x, mean_forces_y, var_disp_x, var_forces_x, var_forces_y);
			assert(save_stats_ok == 0);
		}	
	}

	~SCSimulation_normal() {}
};

#endif /* SCSIMNORMAL_HPP_ */
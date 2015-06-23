#ifndef SCSIMNORMAL_HPP_
#define SCSIMNORMAL_HPP_

#include "UQsim.hpp"
#include "gauss_hermite_quad.hpp"
#include "helper.hpp"

class SCSimulation_normal : public UQSimulation
{
private:
	int rank;
	int ncoeff;
	int nprocs;
	int quad_degree;
	
	GaussHermiteQuadrature ghq;

	int local_global_mapping(const int local_index, const int& rank) const
	{
		int global_id;
		int local_points;

		local_points = ncoeff/nprocs;
		global_id = local_index + rank*local_points;

		return global_id;
	}

public:
	SCSimulation_normal() {}

	SCSimulation_normal(
		std::string& _nastin_dat, 
		std::string& _solidz_dat, 
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

		rho_f_p1 = _rho_f_p1;
		rho_f_p2 = _rho_f_p2;
		nu_f_p1 = _nu_f_p1;
		nu_f_p2 = _nu_f_p2;
		rho_s_p1 = _rho_s_p1;
		rho_s_p2 = _rho_s_p1;
	} 

	virtual int data_decomp() const
	{
		int local_points = 0, rem = 0;

		local_points = ncoeff/nprocs;
		rem = ncoeff % nprocs;

		if (rank == nprocs - 1)   
		{
			local_points += rem;
		}

		return local_points;
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
		std::string modify_nastin_data;
		std::string get_data;
		std::string get_alya_output;
		std::string get_output;

		int run_ok = 0;
		int modify_nastin_data_ok = 0;
		int get_data_ok = 0;
		int gather_alya_output_ok = 0;
		int save_coeff_ok = 0;

		double temp = 0.0;
		double temp_disp_x = 0.0;
		double temp_force0 = 0.0;
		double temp_force1 = 0.0;
		double disp_x = 0.0;
		double force0 = 0.0;
		double force1 = 0.0;

		std::vector<double> disp_x_all;
		std::vector<double> force0_all;
		std::vector<double> force1_all;
		std::vector<double> get_output_values;

		for (int i = 0; i < quad_degree; ++i)
		{
			temp = sqrt(2)*pre_proc_result[i]*nu_f_p2 + nu_f_p1;
			assert(temp >= 0);
			modify_nastin_data = run_insert_nastin_1d(insert_nastin_exec, nastin_dat, temp);

			modify_nastin_data_ok = system(modify_nastin_data.c_str());
			assert(modify_nastin_data_ok >= 0);

			run_ok = system(run_exec.c_str());
			assert(run_ok >= 0);

			get_alya_output = run_gather_alya_output(gather_alya_output, i + 1);
			gather_alya_output_ok = system(get_alya_output.c_str());
			assert(gather_alya_output_ok >=0);

			get_data = run_gather_data(gather_data_exec_sc, output_data, output_file_sc, i + 1);
			get_data_ok = system(get_data.c_str());
			assert(get_data_ok >= 0);

			get_output = run_get_output(get_output_sc, output_file_sc);
			get_output_values = get_output_data(get_output);

			disp_x = get_output_values[0];
			force0 = get_output_values[1];
			force1 = get_output_values[2];

			disp_x_all.push_back(disp_x);
			force0_all.push_back(force0);
			force1_all.push_back(force1);
		}

		for(int j = 0 ; j < ncoeff ; ++j)
		{
			temp_disp_x = 0.0;
			temp_force0 = 0.0;
			temp_force1 = 0.0;

			for (int i = 0; i < quad_degree; ++i)
			{
				temp_disp_x += disp_x_all[i]*ghq.orthogonal_poly(j, sqrt(2)*pre_proc_result[i]) * pre_proc_result[quad_degree + i]/sqrt(M_PI);
				temp_force0 += force0_all[i]*ghq.orthogonal_poly(j, sqrt(2)*pre_proc_result[i]) * pre_proc_result[quad_degree + i]/sqrt(M_PI);
				temp_force1 += force1_all[i]*ghq.orthogonal_poly(j, sqrt(2)*pre_proc_result[i]) * pre_proc_result[quad_degree + i]/sqrt(M_PI);
			}

			temp_disp_x = temp_disp_x/sqrt(ghq.norm_factor(j));
			temp_force0 = temp_force0/sqrt(ghq.norm_factor(j));
			temp_force1 = temp_force1/sqrt(ghq.norm_factor(j));

			save_coeff_ok = save_coeff(coeff_sc, temp_disp_x, temp_force0, temp_force1);
			assert(save_coeff_ok == 1);		
		}
	}

	virtual void post_processing() const
	{
		std::string get_postproc_stat;
		int get_postproc_stat_ok = 0;

		get_postproc_stat = run_postproc_stat(postproc_stat_exec_sc, coeff_sc, postproc_stat_sc);

		get_postproc_stat_ok = system(get_postproc_stat.c_str());
		assert(get_postproc_stat_ok >= 0);
	}

	~SCSimulation_normal() {}
};

#endif /* SCSIMNORMAL_HPP_ */
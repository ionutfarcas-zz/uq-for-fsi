#ifndef SCSIMUNIFORM_HPP_
#define SCSIMUNIFORM_HPP_

#include "UQsim.hpp"
#include "gauss_legendre_quad.hpp"
#include "helper.hpp"

class SCSimulation_uniform : public UQSimulation
{
private:
	int ncoeff;
	int quad_degree;

	GaussLegendreQuadrature glq;

public:
	SCSimulation_uniform() 
	{
		ncoeff = 0;
		quad_degree = 0;
	}

	SCSimulation_uniform(
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
		rho_s_p2 = _rho_s_p2;
	} 

	virtual std::vector<double> pre_processing() const
	{
		std::vector<double> nodes, weights, nodes_weights;

		glq.quad_nodes_weights(quad_degree, weights, nodes);

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

		int no_of_datapoints = 0;
		int no_of_timesteps = 0;
		vec2d_double disp_x_all;
		vec2d_double force0_all;
		vec2d_double force1_all;

		for (int i = 0; i < quad_degree; ++i)
		{
			std::vector<double> disp_x;
			std::vector<double> force0;
			std::vector<double> force1;
			
			temp = (nu_f_p2 - nu_f_p1)/2.0*pre_proc_result[i] + (nu_f_p2 + nu_f_p1)/2.0;
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
			get_output_data(get_output, no_of_datapoints, disp_x, force0, force1);

			disp_x_all.push_back(disp_x);
			force0_all.push_back(force0);
			force1_all.push_back(force1);
		}

		no_of_timesteps = no_of_datapoints/quad_degree;

		for(int timestep = 0 ; timestep < no_of_timesteps ; ++timestep)
		{
			for(int j = 0 ; j < ncoeff ; ++j)
			{
				temp_disp_x = 0.0;
				temp_force0 = 0.0;
				temp_force1 = 0.0;

				for(int i = 0 ; i < quad_degree ; ++i)
				{
					temp_disp_x += disp_x_all[i][timestep]*glq.orthogonal_poly(j, 0.5*pre_proc_result[i] + 0.5) * pre_proc_result[quad_degree + i];
					temp_force0 += force0_all[i][timestep]*glq.orthogonal_poly(j, 0.5*pre_proc_result[i] + 0.5) * pre_proc_result[quad_degree + i];
					temp_force1 += force1_all[i][timestep]*glq.orthogonal_poly(j, 0.5*pre_proc_result[i] + 0.5) * pre_proc_result[quad_degree + i];
				}

				temp_disp_x = (nu_f_p2 - nu_f_p1)*temp_disp_x/glq.norm_factor(j);
				temp_force0 = (nu_f_p2 - nu_f_p1)*temp_force0/glq.norm_factor(j);
				temp_force1 = (nu_f_p2 - nu_f_p1)*temp_force1/glq.norm_factor(j);

				save_coeff_ok = save_coeff(coeff_sc, temp_disp_x, temp_force0, temp_force1);
				assert(save_coeff_ok == 1);		
			}
		}
	}

	virtual void post_processing() const
	{
		std::string get_postproc_stat;
		int get_postproc_stat_ok = 0;

		get_postproc_stat = run_postproc_stat(postproc_stat_exec_sc, coeff_sc, postproc_stat_sc, ncoeff);

		get_postproc_stat_ok = system(get_postproc_stat.c_str());
		assert(get_postproc_stat_ok >= 0);
	}

	~SCSimulation_uniform() {}
};

#endif /* SCSIMUNIFORM_HPP_ */
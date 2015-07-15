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

	int dim;
	
	vec2d_int multi_index_dim2;
	GaussLegendreQuadrature glq;

	int factorial(const int& n) const
	{
		int fact = 0;

		if(n == 0)
		{
			fact = 1;
		}
		else
		{	
			fact = n*factorial(n-1);
		}

		return fact;
	}

	vec2d_int mindex(const int& degree) const
	{
		int j = 0;
		int norm = 0;

		std::vector<int> temp(dim, 0);
		vec2d_int mindex_degree;

		while(true)
		{
			norm = l1_norm(temp);

			if(norm == degree)
			{
				mindex_degree.push_back(temp);
			}

			for(j = dim - 1 ; j >= 0 ; --j)
			{
				if(++temp[j] <= degree)
					break;
				else
					temp[j] = 0;
			}

			if( j < 0)
				break;
		}

		return mindex_degree;
	}

	vec2d_int multi_index() const
	{
		int size = 0;
		vec2d_int m_index_level;
		vec2d_int result;

		for(int i = 0 ; i < ncoeff; ++i)
		{
			m_index_level = this->mindex(i);
			size = m_index_level.size();

			for(int j = size - 1; j >= 0 ; --j)
			{
				result.push_back(m_index_level[j]);
			}
		}

		return result;
	}

public:
	SCSimulation_uniform() {}

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

		dim = 2;

		multi_index_dim2 = multi_index();

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

	int compute_no_coeff() const
	{
		int no_coeff = 1;

		for(int i = 0 ; i < dim ; ++i)
		{
			no_coeff *= (ncoeff + i);
		}

		no_coeff = no_coeff/this->factorial(dim);

		return no_coeff;
	}

	double multi_orthogonal_poly(const double& val1, const double& val2, const int& index) const
	{
		double multi_ortho_poly = 0.0;
		multi_ortho_poly = glq.orthogonal_poly(multi_index_dim2[index][0], val1)*glq.orthogonal_poly(multi_index_dim2[index][1], val2)/(glq.norm_factor(multi_index_dim2[index][0])*glq.norm_factor(multi_index_dim2[index][1]));			

		return multi_ortho_poly;
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
		int ncoeff_2dim = 0;

		std::string modify_nastin_data;
		std::string modify_solidz_data;
		std::string get_data;
		std::string get_alya_output;
		std::string get_output;

		int run_ok = 0;
		int modify_nastin_data_ok = 0;
		int modify_solidz_data_ok = 0;
		int get_data_ok = 0;
		int gather_alya_output_ok = 0;
		int save_coeff_ok = 0;

		double temp_nu_f = 0.0;
		double temp_rho_s = 0.0;
		double temp_disp_x = 0.0;
		double temp_force0 = 0.0;
		double temp_force1 = 0.0;
		double disp_x = 0.0;
		double force0 = 0.0;
		double force1 = 0.0;

		std::vector<double> disp_x_all;
		std::vector<double> force0_all;
		std::vector<double> force1_all;
		vec2d_double disp_x_all_2d;
		vec2d_double force0_all_2d;
		vec2d_double force1_all_2d;

		std::vector<double> get_output_values;
		ncoeff_2dim = this->compute_no_coeff();

		for (int i = 0; i < quad_degree; ++i)
		{
			temp_nu_f = (0.5*pre_proc_result[i] + 0.5)*nu_f_p2 + nu_f_p1;
			assert(temp_nu_f >= 0);
			modify_nastin_data = run_insert_nastin_1d(insert_nastin_exec, nastin_dat, temp_nu_f);
			
			for(int j = 0 ; j < quad_degree ; ++j)
			{
				temp_rho_s = (0.5*pre_proc_result[i] + 0.5)*rho_s_p2 + rho_s_p1;
				assert(temp_rho_s >= 0);
				modify_solidz_data = run_insert_solidz_1d(insert_solidz_exec, solidz_dat, temp_rho_s);

				modify_nastin_data_ok = system(modify_nastin_data.c_str());
				assert(modify_nastin_data_ok >= 0);
				modify_solidz_data_ok = system(modify_solidz_data.c_str());
				assert(modify_solidz_data_ok >= 0);

				run_ok = system(run_exec.c_str());
				assert(run_ok >= 0);
				
				get_alya_output = run_gather_alya_output(gather_alya_output, i*quad_degree + j + 1);
				gather_alya_output_ok = system(get_alya_output.c_str());
				assert(gather_alya_output_ok >=0);

				get_data = run_gather_data(gather_data_exec_sc, output_data, output_file_sc, i*quad_degree + j + 1);
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

			disp_x_all_2d.push_back(disp_x_all);
			force0_all_2d.push_back(force0_all);
			force1_all_2d.push_back(force1_all);
		}

		for(int j = 0 ; j < ncoeff_2dim ; ++j)
		{
			temp_disp_x = 0.0;
			temp_force0 = 0.0;
			temp_force1 = 0.0;

			for (int i = 0; i < quad_degree; ++i)
			{
				for(int k = 0 ; k < quad_degree ; ++k)
				{
					temp_disp_x += disp_x_all_2d[i][k] * multi_orthogonal_poly(0.5*pre_proc_result[i] + 0.5, 0.5*pre_proc_result[k] + 0.5, j) * pre_proc_result[quad_degree + i] * pre_proc_result[quad_degree + k];
					temp_force0 += force0_all_2d[i][k] * multi_orthogonal_poly(0.5*pre_proc_result[i] + 0.5, 0.5*pre_proc_result[k] + 0.5, j) * pre_proc_result[quad_degree + i] * pre_proc_result[quad_degree + k];
					temp_force1 += force1_all_2d[i][k] * multi_orthogonal_poly(0.5*pre_proc_result[i] + 0.5, 0.5*pre_proc_result[k] + 0.5, j) * pre_proc_result[quad_degree + i] * pre_proc_result[quad_degree + k];
				}
			}

			save_coeff_ok = save_coeff(coeff_sc, temp_disp_x, temp_force0, temp_force1);
			assert(save_coeff_ok == 1);		
		}
	}

	virtual void simulation(std::vector<double>& pre_proc_result_dim1, std::vector<double>& pre_proc_result_dim2) const
	{	
		/*TO DO : nothing here */
	}

	virtual void post_processing() const
	{
		std::string get_postproc_stat;
		int get_postproc_stat_ok = 0;

		get_postproc_stat = run_postproc_stat(postproc_stat_exec_sc, coeff_sc, postproc_stat_sc);

		get_postproc_stat_ok = system(get_postproc_stat.c_str());
		assert(get_postproc_stat_ok >= 0);
	}

	~SCSimulation_uniform() {}
};

#endif /* SCSIMUNIFORM_HPP_ */
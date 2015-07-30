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
	int local_quad_points;
	int ncoeff_2dim;

	int dim;
	std::string create_data_each_rank;
	int run_create_data_rank_ok;
	
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

public:
	SCSimulation_uniform() 
	{
		ncoeff = 0;
		quad_degree = 0;
		local_quad_points = 0;
		ncoeff_2dim = 0;

		create_data_each_rank = " ";
		run_create_data_rank_ok = 0;

		dim = 0;
		multi_index_dim2 = {{0}};
	}

	SCSimulation_uniform(std::string& _nastin_dat, 
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

		dim = 2;

		multi_index_dim2 = multi_index();
		ncoeff_2dim = this->compute_no_coeff();

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
		int global_id = 0;

		std::string modify_nastin_data;
		std::string modify_solidz_data;
		int modify_nastin_data_ok = 0;
		int modify_solidz_data_ok = 0;
		
		double temp_nu_f = 0.0;
		double temp_rho_s = 0.0;

		for (int i = 0; i < local_quad_points; ++i)
		{
			global_id = this->local_global_mapping(i, rank);
			
			temp_nu_f = (0.5*pre_proc_result[global_id] + 0.5)*nu_f_p2 + nu_f_p1;
			assert(temp_nu_f >= 0);
			modify_nastin_data = run_insert_nastin_1d(insert_nastin_exec, nastin_dat, temp_nu_f, rank);
			modify_nastin_data_ok = system(modify_nastin_data.c_str());
			assert(modify_nastin_data_ok >= 0);
			
			for(int j = 0 ; j < quad_degree ; ++j)
			{
				temp_rho_s = (0.5*pre_proc_result[j] + 0.5)*rho_s_p2 + rho_s_p1;
				assert(temp_rho_s >= 0);
				modify_solidz_data = run_insert_solidz_1d(insert_solidz_exec, solidz_dat, temp_rho_s, rank);
				modify_solidz_data_ok = system(modify_solidz_data.c_str());
				assert(modify_solidz_data_ok >= 0);
			}
		}
	}

	virtual void simulation(std::vector<double>& pre_proc_result_dim1, std::vector<double>& pre_proc_result_dim2) const
	{	
		/*TO DO : nothing here */
	}

	virtual void post_processing() const
	{
	
	}

	~SCSimulation_uniform() {}
};

#endif /* SCSIMUNIFORM_HPP_ */
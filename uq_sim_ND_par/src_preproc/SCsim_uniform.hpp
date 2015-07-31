#ifndef SCSIMUNIFORM_HPP_
#define SCSIMUNIFORM_HPP_

#include "UQsim.hpp"
#include "gauss_legendre_quad.hpp"

class SCSimulation_uniform : public UQSimulation
{
private:
	int i_ncoeff;
	int i_dim;
	int i_sg_level;
	
	int multi_ncoeff;

	size_t grid_storage_size;	

	SGPP::base::Grid* grid;
	SGPP::base::GridStorage* grid_storage;
	SGPP::base::GridGenerator* grid_gen;
	SGPP::base::OperationQuadrature* quad;

	std::vector<double> l_limits;
	std::vector<double> r_limits;

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

		std::vector<int> temp(i_dim, 0);
		vec2d_int mindex_degree;

		while(true)
		{
			norm = l1_norm(temp);

			if(norm == degree)
			{
				mindex_degree.push_back(temp);
			}

			for(j = i_dim - 1 ; j >= 0 ; --j)
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

		for(int i = 0 ; i < i_ncoeff; ++i)
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

		for(int i = 0 ; i < i_dim ; ++i)
		{
			no_coeff *= (i_ncoeff + i);
		}

		no_coeff = no_coeff/this->factorial(i_dim);

		return no_coeff;
	}

	SGPP::float_t multi_orthogonal_poly(const std::vector<SGPP::float_t>& x, const int& index) const
	{	
		SGPP::float_t multi_ortho_poly = 1.0;
		vec2d_int m_index = this->multi_index();

		for(int j = 0 ; j < i_dim ; ++j)
		{
			multi_ortho_poly *= glq.orthogonal_poly(m_index[index][j], x[j])/glq.norm_factor(m_index[index][j]);
		}

		return multi_ortho_poly;
	}

public:
	SCSimulation_uniform() 
	{
		i_ncoeff = 0;
		i_dim = 0;
		i_sg_level = 0;

		grid_storage_size = 0;
		multi_ncoeff = 0;
		
		nastin_dat = "";
		solidz_dat = "";
		run_exec = "";
		output_data = "";
		gather_data_exec_sc = "";
		get_output_sc = "";
		postproc_stat_exec_sc = "";
		output_file_sc = "";
		coeff_sc = "";
		postproc_stat_sc = "";
		insert_nastin_exec = "";
		insert_solidz_exec = "";
		gather_alya_output = "";

		rho_f_p1 = 0.0;
		rho_f_p2 = 0.0;
		nu_f_p1 = 0.0;
		nu_f_p2 = 0.0;
		rho_s_p1 = 0.0;
		rho_s_p2 = 0.0;

		grid = nullptr;
		grid_storage = nullptr;
		grid_gen = nullptr;
		quad = nullptr;
	}

	SCSimulation_uniform(
		std::string& _nastin_dat, 
		std::string& _solidz_dat,
		std::string& _create_data_point,  
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
		const unsigned int& _prob_dim,
		const unsigned int& _sg_level,
		const double& _rho_f_p1, 
		const double& _rho_f_p2, 
		const double& _nu_f_p1, 
		const double& _nu_f_p2, 
		const double& _rho_s_p1, 
		const double& _rho_s_p2)
	{
		i_ncoeff = _ncoeff;
		i_dim = _prob_dim;
		i_sg_level = _sg_level;

		nastin_dat = _nastin_dat;
		solidz_dat = _solidz_dat;
		create_data_point = _create_data_point;
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

		l_limits = {_nu_f_p1, _rho_s_p1};
		r_limits = {_nu_f_p2, _rho_s_p2};

		grid = SGPP::base::Grid::createLinearGrid(_prob_dim);
		grid_storage = grid->getStorage(); 
		grid_gen = grid->createGridGenerator();
		grid_gen->regular(_sg_level);

		quad = SGPP::op_factory::createOperationQuadrature(*grid);

		grid_storage_size = grid_storage->size();
		multi_ncoeff = compute_no_coeff();
	}

	size_t get_storage_size() const
	{
		return this->grid_storage_size;
	}

	virtual double compute_volume() const
	{
		double volume = 1.0;

		for(int i = 0 ; i < i_dim ; ++i)
		{
			volume *= (r_limits[i] - l_limits[i]);
		}

		return volume;
	}

	virtual vec2d_float_t pre_processing() const
	{
		SGPP::base::GridIndex* gp;
		double grid_point = 0.0;
		vec2d_float_t result;

		for(size_t i = 0; i < grid_storage_size; ++i) 
		{
			std::vector<SGPP::float_t> p;
			gp = grid_storage->get(i);
			
			for(int j = 0 ; j < i_dim ; ++j)
			{
				grid_point = l_limits[j] + (r_limits[j] - l_limits[j])*gp->getCoord(j);
				p.push_back(grid_point);
			}

			result.push_back(p);
		}

		return result;
	}

	virtual std::vector<double> pre_processing(const double& param1, const double& param2) const
	{
		// do nothing; useful only for MCS
		std::vector<double> dummy;

		return dummy;
	}

	virtual void simulation(const vec2d_float_t& pre_proc_result) const
	{
		std::string modify_nastin_data;
		std::string modify_solidz_data;
		int modify_nastin_data_ok = 0;
		int modify_solidz_data_ok = 0;
		
		double temp_nu_f = 0.0;
		double temp_rho_s = 0.0;

		std::string create_data_each_point;
		int run_create_data_point_ok = 0;

		for(size_t i = 0; i < grid_storage_size; ++i) 
		{	
			create_data_each_point = run_create_data_point(create_data_point, i);
			run_create_data_point_ok = system(create_data_each_point.c_str());
			assert(run_create_data_point_ok >=0 );

			temp_nu_f = (0.5*pre_proc_result[i][0] + 0.5)*nu_f_p2 + nu_f_p1;
			assert(temp_nu_f >= 0);
			modify_nastin_data = run_insert_nastin_1d(insert_nastin_exec, nastin_dat, temp_nu_f, i);

			temp_rho_s = (0.5*pre_proc_result[i][0] + 0.5)*rho_s_p2 + rho_s_p1;
			assert(temp_rho_s >= 0);
			modify_solidz_data = run_insert_solidz_1d(insert_solidz_exec, solidz_dat, temp_rho_s, i);

			modify_nastin_data_ok = system(modify_nastin_data.c_str());
			assert(modify_nastin_data_ok >= 0);
			modify_solidz_data_ok = system(modify_solidz_data.c_str());
			assert(modify_solidz_data_ok >= 0);
		}
	}

	virtual void simulation(std::vector<double>& pre_proc_result_dim1, std::vector<double>& pre_proc_result_dim2) const
	{
		/*TO DO - Nothing here */
	}

	virtual void post_processing() const
	{
		
	}

	~SCSimulation_uniform() {}
};

#endif /* SCSIMUNIFORM_HPP_ */
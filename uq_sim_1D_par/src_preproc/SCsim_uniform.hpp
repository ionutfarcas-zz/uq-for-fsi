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

	std::string create_data_each_rank;
	int run_create_data_rank_ok;

	GaussLegendreQuadrature glq;

public:
	SCSimulation_uniform() 
	{
		ncoeff = 0;
		quad_degree = 0;
		local_quad_points = 0;

		create_data_each_rank = " ";
		run_create_data_rank_ok = 0;
	}

	SCSimulation_uniform(
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
		int modify_nastin_data_ok = 0;
		
		double temp = 0.0;
		
		for (int i = 0; i < local_quad_points; ++i)
		{
			global_id = this->local_global_mapping(i, rank);
			
			temp = (nu_f_p2 - nu_f_p1)/2.0*pre_proc_result[global_id] + (nu_f_p2 + nu_f_p1)/2.0;
			assert(temp >= 0);
			modify_nastin_data = run_insert_nastin_1d(insert_nastin_exec, nastin_dat, temp, rank);
			modify_nastin_data_ok = system(modify_nastin_data.c_str());
			assert(modify_nastin_data_ok >= 0);
		}
	}

	virtual void post_processing() const
	{
		
	}

	~SCSimulation_uniform() {}
};

#endif /* SCSIMUNIFORM_HPP_ */
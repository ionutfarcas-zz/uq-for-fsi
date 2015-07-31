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
	
	GaussHermiteQuadrature ghq;

public:
	SCSimulation_normal() 
	{
		ncoeff = 0;
		quad_degree = 0;
	}

	SCSimulation_normal(
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
		int modify_nastin_data_ok = 0;
		
		double temp = 0.0;

		std::string create_data_each_point;
		int run_create_data_point_ok = 0;

		for (int i = 0; i < quad_degree; ++i)
		{
			create_data_each_point = run_create_data_point(create_data_point, i);
			run_create_data_point_ok = system(create_data_each_point.c_str());
			assert(run_create_data_point_ok >=0 );

			temp = sqrt(2.0)*pre_proc_result[i]*nu_f_p2 + nu_f_p1;
			assert(temp >= 0);
			modify_nastin_data = run_insert_nastin_1d(insert_nastin_exec, nastin_dat, temp, i);
			modify_nastin_data_ok = system(modify_nastin_data.c_str());
			assert(modify_nastin_data_ok >= 0);
		}
	}

	virtual void post_processing() const
	{
		
	}

	~SCSimulation_normal() {}
};

#endif /* SCSIMNORMAL_HPP_ */
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
	int local_points;	

	SGPP::base::Grid* grid;
	SGPP::base::GridStorage* grid_storage;
	SGPP::base::GridGenerator* grid_gen;
	SGPP::base::OperationQuadrature* quad;

	std::vector<double> l_limits;
	std::vector<double> r_limits;

	std::string create_data_each_rank;
	int run_create_data_rank_ok;

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

		local_points = 0;

		create_data_each_rank = " ";
		run_create_data_rank_ok = 0;

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
		const unsigned int& _prob_dim,
		const unsigned int& _sg_level,
		unsigned int _rank,
		const unsigned int& _nprocs,
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

		rank = _rank;
		nprocs = _nprocs;

		rho_f_p1 = _rho_f_p1;
		rho_f_p2 = _rho_f_p2;
		nu_f_p1 = _nu_f_p1;
		nu_f_p2 = _nu_f_p2;
		rho_s_p1 = _rho_s_p1;
		rho_s_p2 = _rho_s_p2;

		l_limits = {_nu_f_p1, _rho_s_p1};
		r_limits = {_nu_f_p2, _rho_s_p2};

		create_data_each_rank = run_create_data_rank(create_data_rank, _rank);
		run_create_data_rank_ok = system(create_data_each_rank.c_str());
		assert(run_create_data_rank_ok >=0 );

		grid = SGPP::base::Grid::createLinearGrid(_prob_dim);
		grid_storage = grid->getStorage(); 
		grid_gen = grid->createGridGenerator();
		grid_gen->regular(_sg_level);

		quad = SGPP::op_factory::createOperationQuadrature(*grid);

		grid_storage_size = grid_storage->size();
		local_points = data_decomp();
		multi_ncoeff = compute_no_coeff();
	}

	size_t get_storage_size() const
	{
		return this->grid_storage_size;
	}

	virtual int local_global_mapping(const int local_index, const int& rank) const
	{
		int global_id;
		int local_points_proc;

		local_points_proc = grid_storage_size/nprocs;
		global_id = local_index + rank*local_points_proc;

		return global_id;
	}

	virtual int data_decomp() const
	{
		int local_points_proc = 0;
		int rem = 0;

		local_points_proc = grid_storage_size/nprocs;
		rem = grid_storage_size % nprocs;

		if (rank == nprocs - 1)   
		{
			local_points_proc += rem;
		}

		return local_points_proc;
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
		int global_id = 0;
		double volume = 0.0;

		std::string modify_nastin_data;
		std::string modify_solidz_data;
		std::string alya_nastin_solidz;
		std::string get_data;
		std::string get_alya_output;
		std::string get_output;

		int run_alya_ok = 0;
		int modify_nastin_data_ok = 0;
		int modify_solidz_data_ok = 0;
		int get_data_ok = 0;
		int gather_alya_output_ok = 0;

		int no_valid_lines_local = 0;
		int no_valid_lines = 0;

		double temp_nu_f = 0.0;
		double temp_rho_s = 0.0;
		double temp_disp_x = 0.0;
		double temp_force0 = 0.0;
		double temp_force1 = 0.0;
		double disp_x_local = 0.0;
		double force0_local = 0.0;
		double force1_local = 0.0;

		vec2d_double get_output_values;

		int* scount = (int*) malloc(nprocs * sizeof(int));
		int* displs = (int*) malloc(nprocs * sizeof(int));

		double* alpha_xdisp_local = (double*)calloc(local_points, sizeof(double));
		double* alpha_force0_local = (double*)calloc(local_points, sizeof(double));
		double* alpha_force1_local = (double*)calloc(local_points, sizeof(double));

		double* alpha_xdisp_temp;
		double* alpha_force0_temp;
		double* alpha_force1_temp;

		if(rank == 0)
		{
			alpha_xdisp_temp = (double*)calloc(grid_storage_size, sizeof(double));
			alpha_force0_temp = (double*)calloc(grid_storage_size, sizeof(double));
			alpha_force1_temp = (double*)calloc(grid_storage_size, sizeof(double));
		}

		for (int i = 0 ; i < nprocs; ++i)
		{
			scount[i] = local_points;
			displs[i] = i * local_points;
		}

		volume = this->compute_volume();

		for(int i = 0; i < local_points; ++i) 
		{	
			global_id = this->local_global_mapping(i, rank);	

			temp_nu_f =  sqrt(2)*pre_proc_result[global_id][0]*nu_f_p2 + nu_f_p1;
			assert(temp_nu_f >= 0);
			modify_nastin_data = run_insert_nastin_1d(insert_nastin_exec, nastin_dat, temp_nu_f, rank);

			temp_rho_s = sqrt(2)*pre_proc_result[global_id][1]*rho_s_p2 + rho_s_p1;
			assert(temp_rho_s >= 0);
			modify_solidz_data = run_insert_solidz_1d(insert_solidz_exec, solidz_dat, temp_rho_s, rank);

			modify_nastin_data_ok = system(modify_nastin_data.c_str());
			assert(modify_nastin_data_ok >= 0);
			modify_solidz_data_ok = system(modify_solidz_data.c_str());
			assert(modify_solidz_data_ok >= 0);

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

			alpha_xdisp_local[i] = disp_x_local;
			alpha_force0_local[i] = force0_local;
			alpha_force1_local[i] = force1_local;
		}

		MPI_Gatherv(alpha_xdisp_local, local_points, MPI_DOUBLE, alpha_xdisp_temp, scount, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gatherv(alpha_force0_local, local_points, MPI_DOUBLE, alpha_force0_temp, scount, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gatherv(alpha_force1_local, local_points, MPI_DOUBLE, alpha_force1_temp, scount, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		MPI_Allreduce(&no_valid_lines_local, &no_valid_lines, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		/* Insert fault-tolerant mechanism here by reducing the quadrature degre */

		if(rank == 0)
		{	
			int save_coeff_ok = 0;

			SGPP::base::DataVector alpha_xdisp(grid_storage_size);
			SGPP::base::DataVector alpha_force0(grid_storage_size);
			SGPP::base::DataVector alpha_force1(grid_storage_size);

			for(int k = 0 ; k < multi_ncoeff ; ++k)
			{	
				alpha_xdisp.setAll(0.0);
				alpha_force0.setAll(0.0);
				alpha_force1.setAll(0.0);

				for(size_t i = 0; i < grid_storage_size; ++i) 
				{
					alpha_xdisp[i] = alpha_xdisp_temp[i]*multi_orthogonal_poly(pre_proc_result[i], k);
					alpha_force0[i] = alpha_force0_temp[i]*multi_orthogonal_poly(pre_proc_result[i], k);
					alpha_force1[i] = alpha_force1_temp[i]*multi_orthogonal_poly(pre_proc_result[i], k);
				}
				
				SGPP::op_factory::createOperationHierarchisation(*grid)->doHierarchisation(alpha_xdisp);
				SGPP::op_factory::createOperationHierarchisation(*grid)->doHierarchisation(alpha_force0);
				SGPP::op_factory::createOperationHierarchisation(*grid)->doHierarchisation(alpha_force1);

				temp_disp_x = volume*quad->doQuadrature(alpha_xdisp);
				temp_force0 = volume*quad->doQuadrature(alpha_force0);
				temp_force1 = volume*quad->doQuadrature(alpha_force1);

				save_coeff_ok = save_coeff(coeff_sc, temp_disp_x, temp_force0, temp_force1);
				assert(save_coeff_ok == 1);	
			}

			free(alpha_xdisp_temp);
			free(alpha_force0_temp);
			free(alpha_force1_temp);
		}

		free(scount);
		free(displs);
		free(alpha_xdisp_local);
		free(alpha_force0_local);
		free(alpha_force1_local);
	}

	virtual void simulation(std::vector<double>& pre_proc_result_dim1, std::vector<double>& pre_proc_result_dim2) const
	{
		/*TO DO - Nothing here */
	}

	virtual void post_processing() const
	{
		if(rank == 0)
		{
			int get_coeff_ok = 0;
			int save_stats_ok = 0;

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

			for(int i = 1 ; i < multi_ncoeff ; ++i)
			{
				var_disp_x += disp_x[i]*disp_x[i];
				var_forces_x += force0[i]*force0[i];
				var_forces_y += force1[i]*force1[i];
			}

			save_stats_ok = write_stat_to_file(postproc_stat_sc, mean_disp_x, mean_forces_x, mean_forces_y, var_disp_x, var_forces_x, var_forces_y);
			assert(save_stats_ok == 0);
		}	
	}

	~SCSimulation_uniform() {}
};

#endif /* SCSIMUNIFORM_HPP_ */
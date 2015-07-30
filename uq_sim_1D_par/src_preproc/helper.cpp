#include "helper.hpp"

template<typename T>
T str_to_number(const std::string& no)
{
	T value;
	std::stringstream stream(no);
	stream >> value;

	if (stream.fail()) 
	{
		std::runtime_error e(no);
		throw e;
	}

	return value;
}

std::string run_insert_nastin_1d(const std::string nastin_exec, std::string nastin_data, 
	double& new_viscosity, int rank)
{
	std::string delimiter = "/";
	std::string token;
	std::string new_nastin_data;
	std::string insert_rank = "_" + std::to_string(rank);
	size_t pos = 0;
	int iter = 0;

	std::stringstream caller;

	while ((pos = nastin_data.find(delimiter)) != std::string::npos) 
	{
		token = nastin_data.substr(0, pos);
		new_nastin_data += token;
		
		if(iter == 2 || iter == 3)
		{
			new_nastin_data += insert_rank;
		}

		new_nastin_data += delimiter;
		nastin_data.erase(0, pos + delimiter.length());
		++iter;
	}

	new_nastin_data += nastin_data;
	caller << nastin_exec << " " << new_nastin_data << " " << new_viscosity;

	return caller.str();
}

std::string run_insert_nastin_2d(const std::string nastin_exec, std::string nastin_data, 
	double& new_density, double& new_viscosity, int rank)
{
	std::string delimiter = "/";
	std::string token;
	std::string new_nastin_data;
	std::string insert_rank = "_" + std::to_string(rank);
	size_t pos = 0;
	int iter = 0;

	std::stringstream caller;

	while ((pos = nastin_data.find(delimiter)) != std::string::npos) 
	{
		token = nastin_data.substr(0, pos);
		new_nastin_data += token;
		
		if(iter == 2 || iter == 3)
		{
			new_nastin_data += insert_rank;
		}

		new_nastin_data += delimiter;
		nastin_data.erase(0, pos + delimiter.length());
		++iter;
	}

	new_nastin_data += nastin_data;
	caller << nastin_exec << " " << new_nastin_data << " " << new_density << " " << new_viscosity;

	return caller.str();
}

std::string run_insert_solidz_1d(const std::string solidz_exec, std::string solidz_data, 
	double& new_density, int rank)
{
	std::string delimiter = "/";
	std::string token;
	std::string new_solidz_data;
	std::string insert_rank = "_" + std::to_string(rank);
	size_t pos = 0;
	int iter = 0;

	std::stringstream caller;

	while ((pos = solidz_data.find(delimiter)) != std::string::npos) 
	{
		token = solidz_data.substr(0, pos);
		new_solidz_data += token;
		
		if(iter == 2 || iter == 3)
		{
			new_solidz_data += insert_rank;
		}

		new_solidz_data += delimiter;
		solidz_data.erase(0, pos + delimiter.length());
		++iter;
	}

	new_solidz_data += solidz_data;
	caller << solidz_exec << " " << new_solidz_data << " " << new_density;

	return caller.str();
}

std::string run_alya(const std::string& alya_run, const int& rank)
{
	std::stringstream caller;
	caller << alya_run << " " << rank;

	return caller.str();
}

std::string run_create_data_rank(const std::string create_data_rank, const int& rank)
{
	std::stringstream caller;
	caller << create_data_rank << " " << rank;

	return caller.str();
}

std::string run_gather_data(const std::string gather_data_exec, std::string datafile, std::string datafile_all, const int& id, int rank)
{
	std::string delimiter = "/";
	std::string token;
	std::string new_datafile;
	std::string new_datafile_all;
	std::string insert_rank = "_" + std::to_string(rank);
	std::string rank_insert = std::to_string(rank) + "_";
	size_t pos = 0;
	int iter_data = 0;

	std::stringstream caller;

	while ((pos = datafile.find(delimiter)) != std::string::npos) 
	{
		token = datafile.substr(0, pos);
		new_datafile += token;
		
		if(iter_data == 2 || iter_data == 3)
		{
			new_datafile += insert_rank;
		}

		new_datafile += delimiter;
		datafile.erase(0, pos + delimiter.length());
		++iter_data;
	}
	new_datafile += datafile;

	while ((pos = datafile_all.find(delimiter)) != std::string::npos) 
	{
		token = datafile_all.substr(0, pos);
		new_datafile_all += token;

		new_datafile_all += delimiter;
		datafile_all.erase(0, pos + delimiter.length());
	}
	new_datafile_all += rank_insert + datafile_all;

	caller << gather_data_exec << " " << new_datafile << " " << new_datafile_all << " " << id;

	return caller.str();
}

std::string run_postproc_stat(const std::string postproc_stat, std::string datafile_all, const int& no_of_simulations, int rank)
{
	std::string delimiter = "/";
	std::string token;
	std::string new_datafile_all;
	std::string rank_insert = std::to_string(rank) + "_";
	size_t pos = 0;

	std::stringstream caller;

	while ((pos = datafile_all.find(delimiter)) != std::string::npos) 
	{
		token = datafile_all.substr(0, pos);
		new_datafile_all += token;

		new_datafile_all += delimiter;
		datafile_all.erase(0, pos + delimiter.length());
	}
	new_datafile_all += rank_insert + datafile_all;

	caller << postproc_stat << " " << new_datafile_all << " " << no_of_simulations;

	return caller.str();
}

std::vector<double> get_stat_mc(const std::string& get_sums, int& no_of_datapoints)
{
	FILE* stream;
	char buffer[512];

	std::string no_of_datapoints_str;
	std::string sum_disp_x_str;
	std::string sum_forces_x_str;
	std::string sum_forces_y_str;
	std::string sum2_disp_x_str;
	std::string sum2_forces_x_str;
	std::string sum2_forces_y_str;

	double sum_disp_x = 0.0;
	double sum_forces_x = 0.0;
	double sum_forces_y = 0.0;
	double sum2_disp_x = 0.0;
	double sum2_forces_x = 0.0;
	double sum2_forces_y = 0.0;

	std::vector<double> results;

	stream = popen(get_sums.c_str(), "r");

	if(stream) 
	{
		while (!feof(stream))
		{
			if (fgets(buffer, sizeof(buffer), stream) != NULL)
			{
				std::stringstream temp(buffer);
				temp >> no_of_datapoints_str >> sum_disp_x_str >> sum_forces_x_str >> sum_forces_y_str >> sum2_disp_x_str >> sum2_forces_x_str >> sum2_forces_y_str;

				no_of_datapoints = str_to_number<int>(no_of_datapoints_str);

				sum_disp_x = str_to_number<double>(sum_disp_x_str);
				sum_forces_x = str_to_number<double>(sum_forces_x_str);
				sum_forces_y = str_to_number<double>(sum_forces_y_str);
				sum2_disp_x = str_to_number<double>(sum2_disp_x_str);
				sum2_forces_x = str_to_number<double>(sum2_forces_x_str);
				sum2_forces_y = str_to_number<double>(sum2_forces_y_str);

				results.push_back(sum_disp_x);
				results.push_back(sum_forces_x);
				results.push_back(sum_forces_y);
				results.push_back(sum2_disp_x);
				results.push_back(sum2_forces_x);
				results.push_back(sum2_forces_y);
			}
		}

		pclose(stream);
	}
	else
	{
		throw "Error reading script output!"; 
	}

	return results;
}

int write_stat_to_file(
	const std::string& stat_file_name, 
	const double& mean_disp_x, 
	const double& mean_forces_x, 
	const double& mean_forces_y, 
	const double& var_disp_x, 
	const double& var_forces_x, 
	const double& var_forces_y,
	const int& timestep)
{
	std::ofstream stat_file;

	stat_file.open(stat_file_name.c_str(), std::ios::app);

	if (!stat_file.is_open())
	{
		std::cout << "Could not open config file!" << std::endl;
		return 0;
	}

	stat_file << "Timestep: " << timestep << " The mean of the displacement on the x axis is: " << mean_disp_x << std::endl;
	stat_file << "Timestep: " << timestep << " The variance of the displacement on the x axis is: " << var_disp_x << std::endl;

	stat_file << "Timestep: " << timestep << " The mean of force x is: " << mean_forces_x << std::endl;
	stat_file << "Timestep: " << timestep << " The variance of the force x  is: " << var_forces_x << std::endl;
	
	stat_file << "Timestep: " << timestep << " The mean of force y is: " << mean_forces_y << std::endl;
	stat_file << "Timestep: " << timestep << " The variance of the force y  is: " << var_forces_y << std::endl;
	
	stat_file << std::endl;

	stat_file.close();

  	return 0;
}

std::string run_get_output(const std::string get_output, std::string data, int rank)
{
	std::string delimiter = "/";
	std::string token;
	std::string new_data;
	std::string rank_insert = std::to_string(rank) + "_";
	size_t pos = 0;

	std::stringstream caller;

	while ((pos = data.find(delimiter)) != std::string::npos) 
	{
		token = data.substr(0, pos);
		new_data += token;

		new_data += delimiter;
		data.erase(0, pos + delimiter.length());
	}
	new_data += rank_insert + data;

	caller << get_output << " " << new_data;

	return caller.str();
}

int get_coeff_sc(const std::string& coeff_sc, 
	std::vector<double>& disp_x, 
	std::vector<double>& force0, 
	std::vector<double>& force1)
{
	std::ifstream coeff_file;

	std::string line;
	double disp_x_val = 0.0;
	double force0_val = 0.0;
	double force1_val = 0.0;
	std::string disp_x_val_str;
	std::string forces0_val_str;
	std::string forces1_val_str;

	coeff_file.open(coeff_sc.c_str(), std::ios::in);

	if (!coeff_file.is_open())
	{
		std::cout << "Could not open coeff_file file!" << std::endl;
		return 0;
	}

	while (std::getline(coeff_file, line))
    {
        std::stringstream oneline(line);

        oneline >> disp_x_val_str >> forces0_val_str >> forces1_val_str;
        disp_x_val = str_to_number<double>(disp_x_val_str);
        force0_val = str_to_number<double>(forces0_val_str);
        force1_val = str_to_number<double>(forces1_val_str);

        disp_x.push_back(disp_x_val);
        force0.push_back(force0_val);
        force1.push_back(force1_val);
    }

    coeff_file.close();

    return 1;
}

std::string run_gather_alya_output(const std::string get_alya_output, const int& run_id)
{
	std::stringstream caller;
	caller << get_alya_output << " " << run_id;

	return caller.str();
}

int parse_configfile(const std::string& config_file_name,
	std::string& nastin_dat,
	std::string& solidz_dat,
	std::string& create_data_rank,
	std::string& run_exec,
	std::string& output_data,
	std::string& gather_data_exec_mc,
	std::string& postproc_stat_exec_mc,
	std::string& postproc_file_all_mc,
	std::string& postproc_stat_mc,
	std::string& gather_data_exec_sc,
	std::string& get_output_sc,
	std::string& postproc_stat_exec_sc,
	std::string& output_file_sc,
	std::string& coeff_sc,
	std::string& postproc_stat_sc,
	std::string& insert_nastin_exec,
	std::string& insert_solidz_exec,
	std::string& gather_alya_output,
	unsigned int& uq_method,
	unsigned int& pdf,
	unsigned int& nsamples,
	unsigned int& ncoeff,
	unsigned int& quad_degree,
	double& rho_f_p1,
	double& rho_f_p2,
	double& nu_f_p1,
	double& nu_f_p2,
	double& rho_s_p1,
	double& rho_s_p2)
{
	std::ifstream config_file;
	std::string line, token_1, token_2, token_3;

	config_file.open(config_file_name.c_str(), std::ios::in);

	if (!config_file.is_open())
	{
		std::cout << "Could not open config file!" << std::endl;
		return 0;
	}

	while(getline(config_file, line))
	{
		if(line.at(0) != '#')
		{
			std::stringstream tokenizer(line);

			tokenizer >> token_1 >> token_2 >> token_3;

			if(token_1.compare("nastin_dat") == 0)
			{
				nastin_dat = token_3;
			}
			else if(token_1.compare("solidz_dat") == 0)
			{
				solidz_dat = token_3;
			}
			else if(token_1.compare("create_data_rank") == 0)
			{
				create_data_rank = token_3;
			}
			else if(token_1.compare("run_exec") == 0)
			{
				run_exec = token_3;
			}
			else if(token_1.compare("output_data") == 0)
			{
				output_data = token_3;
			}
			else if(token_1.compare("gather_data_exec_mc") == 0)
			{
				gather_data_exec_mc = token_3;
			}
			else if(token_1.compare("postproc_stat_exec_mc") == 0)
			{
				postproc_stat_exec_mc = token_3;
			}
			else if(token_1.compare("postproc_file_all_mc") == 0)
			{
				postproc_file_all_mc = token_3;
			}
			else if(token_1.compare("postproc_stat_mc") == 0)
			{
				postproc_stat_mc = token_3;
			}
			else if(token_1.compare("gather_data_exec_sc") == 0)
			{
				gather_data_exec_sc = token_3;
			}
			else if(token_1.compare("get_output_sc") == 0)
			{
				get_output_sc = token_3;
			}
			else if(token_1.compare("postproc_stat_exec_sc") == 0)
			{
				postproc_stat_exec_sc = token_3;
			}
			else if(token_1.compare("output_file_sc") == 0)
			{
				output_file_sc = token_3;
			}
			else if(token_1.compare("coeff_sc") == 0)
			{
				coeff_sc = token_3;
			}
			else if(token_1.compare("postproc_stat_sc") == 0)
			{
				postproc_stat_sc = token_3;
			}
			else if(token_1.compare("insert_nastin_exec") == 0)
			{
				insert_nastin_exec = token_3;
			}
			else if(token_1.compare("insert_solidz_exec") == 0)
			{
				insert_solidz_exec = token_3;
			}
			else if(token_1.compare("gather_alya_output") == 0)
			{
				gather_alya_output = token_3;
			}
			else if(token_1.compare("uq_method") == 0)
			{
				uq_method = str_to_number<unsigned int>(token_3);
			}
			else if(token_1.compare("pdf") == 0)
			{
				pdf = str_to_number<unsigned int>(token_3);
			}
			else if(token_1.compare("nsamples") == 0)
			{
				nsamples = str_to_number<unsigned int>(token_3);
			}
			else if(token_1.compare("ncoeff") == 0)
			{
				ncoeff = str_to_number<unsigned int>(token_3);
			}
			else if(token_1.compare("quad_degree") == 0)
			{
				quad_degree = str_to_number<unsigned int>(token_3);
			}
			else if(token_1.compare("rho_f_p1") == 0)
			{
				rho_f_p1 = str_to_number<double>(token_3);
			}
			else if(token_1.compare("rho_f_p2") == 0)
			{
				rho_f_p2 = str_to_number<double>(token_3);
			}
			else if(token_1.compare("nu_f_p1") == 0)
			{
				nu_f_p1 = str_to_number<double>(token_3);
			}
			else if(token_1.compare("nu_f_p2") == 0)
			{
				nu_f_p2 = str_to_number<double>(token_3);
			}
			else if(token_1.compare("rho_s_p1") == 0)
			{
				rho_s_p1 = str_to_number<double>(token_3);
			}
			else if(token_1.compare("rho_s_p2") == 0)
			{
				rho_s_p2 = str_to_number<double>(token_3);
			}
			else
			{
				std::cout << "Error processing config file! Some parameters are missing" << std::endl;
			}
		}
	}

	config_file.close();

	return 1;	
}

vec2d_double get_output_data(const std::string get_output_sc, int& no_of_datapoints)
{
	FILE* stream;
	char buffer[256];
	std::string no_valid_lines_str;
	std::string disp_x_str;
	std::string force0_str;
	std::string force1_str;
	double disp_x;
	double force0;
	double force1;

	vec2d_double data_all;

	stream = popen(get_output_sc.c_str(), "r");
	if (stream) 
	{
		while (!feof(stream))
		{
			if (fgets(buffer, sizeof(buffer), stream) != NULL)
			{
				std::vector<double> data;
				std::stringstream temp(buffer);
				temp >> disp_x_str >> force0_str >> force1_str;

				++no_of_datapoints;
				disp_x = str_to_number<double>(disp_x_str);
				force0 = str_to_number<double>(force0_str);
				force1 = str_to_number<double>(force1_str);

				data.push_back(disp_x);
				data.push_back(force0);
				data.push_back(force1);

				data_all.push_back(data);
			}
		}
		pclose(stream);
	}
	else
	{
		throw "Error reading data file!"; 
	}

	return data_all;
}

int save_coeff(const std::string file_name, const double& disp_x, const double& force0, const double& force1)
{
	std::ofstream coeff_file; 

	coeff_file.open(file_name.c_str(), std::ios::app);

	if (coeff_file.is_open())
	{
		coeff_file << disp_x << " " << force0 << " " << force1 << std::endl;
		coeff_file.close();

		return 1;	
	}
	else
	{
		std::cout << "Could not open coefficiants file!" << std::endl;

		return 0;
	}
}

void Simulation_Message(const std::string& txt)
{
	int myrank;

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	std::cout << "-MESSAGE- P:" << myrank << " " << txt << std::endl;
	fflush(stdout);
	fflush(stderr);
}


void Simulation_Sync(const std::string& txt)
{
	int myrank;

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Barrier(MPI_COMM_WORLD);                               
	std::cout << "-MESSAGE- P:" << myrank << " " << txt << std::endl;
	fflush(stdout);
	fflush(stderr);
	MPI_Barrier(MPI_COMM_WORLD);
}


void Simulation_Stop(const std::string& txt)
{
	int myrank;

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Barrier(MPI_COMM_WORLD);
	std::cout << "-STOP- P:" << myrank << " " << txt << std::endl;
	fflush(stdout);
	fflush(stderr);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}
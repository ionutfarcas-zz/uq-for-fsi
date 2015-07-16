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

std::string run_insert_nastin_1d(const std::string nastin_exec, const std::string nastin_data, 
	double& new_viscosity)
{
	std::stringstream caller;
	caller << nastin_exec << " " << nastin_data << " " << new_viscosity;

	return caller.str();
}

std::string run_insert_nastin_2d(const std::string nastin_exec, const std::string nastin_data, 
	double& new_density, double& new_viscosity)
{
	std::stringstream caller;
	caller << nastin_exec << " " << nastin_data << " " << new_density << " " << new_viscosity;

	return caller.str();
}

std::string run_insert_solidz_1d(const std::string solidz_exec, const std::string solidz_data, 
	double& new_density)
{
	std::stringstream caller;
	caller << solidz_exec << " " << solidz_data << " " << new_density;

	return caller.str();
}

std::string run_gather_data(const std::string gather_data_exec, const std::string datafile, 
	const std::string datafile_all, const int& id)
{
	std::stringstream caller;
	caller << gather_data_exec << " " << datafile << " " << datafile_all << " " << id;

	return caller.str();
}

std::string run_postproc_stat(const std::string postproc_stat, const std::string all_data, 
	const std::string stats)
{
	std::stringstream caller;
	caller << postproc_stat << " " << all_data << " " << stats;

	return caller.str();
}

std::string run_get_output(const std::string get_output, const std::string data)
{
	std::stringstream caller;
	caller << get_output << " " << data;

	return caller.str();
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

std::vector<double> get_output_data(const std::string get_output_sc) 
{
    FILE* stream;
    char buffer[256];
    std::string disp_x_str;
    std::string force0_str;
    std::string force1_str;
    double disp_x;
    double force0;
    double force1;

    std::vector<double> data;

    stream = popen(get_output_sc.c_str(), "r");
    if (stream) 
    {
        while (!feof(stream))
        {
            if (fgets(buffer, sizeof(buffer), stream) != NULL)
            {
                std::stringstream temp(buffer);
                temp >> disp_x_str >> force0_str >> force1_str;
                
                disp_x = str_to_number<double>(disp_x_str);
                force0 = str_to_number<double>(force0_str);
                force1 = str_to_number<double>(force1_str);

                data.push_back(disp_x);
                data.push_back(force0);
                data.push_back(force1);
            }
        }
        pclose(stream);
    }
    else
    {
        throw "Error reading data file!"; 
    }

    return data;
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
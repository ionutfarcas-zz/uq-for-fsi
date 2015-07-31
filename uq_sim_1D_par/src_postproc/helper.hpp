#ifndef HELPER_HPP_
#define HELPER_HPP_

#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <valarray>
#include <chrono>

typedef std::vector<std::vector<double>> vec2d_double;

std::string run_insert_nastin_2d(const std::string nastin_exec, std::string nastin_data, 
	double& new_density, double& new_viscosity, int point);

std::string run_insert_nastin_1d(const std::string nastin_exec, std::string nastin_data, 
	double& new_viscosity, int point);

std::string run_insert_solidz_1d(const std::string solidz_exec, std::string solidz_data, 
	double& new_density, int point);

std::string run_alya(const std::string& alya_run, const int& point);

std::string run_create_data_point(const std::string create_data_point, const int& point);

std::string run_gather_data(const std::string gather_data_exec, std::string datafile, std::string datafile_all, const int& id, int point);

std::string run_postproc_stat(const std::string postproc_stat, std::string datafile_all, const int& no_of_simulations, int point);

std::vector<double> get_stat_mc(const std::string& get_sums, int& no_of_datapoints);

int write_stat_to_file(
	const std::string& stat_file_name, 
	const double& mean_disp_x, 
	const double& mean_forces_x, 
	const double& mean_forces_y, 
	const double& stddev_disp_x, 
	const double& stddev_forces_x, 
	const double& stddev_forces_y,
	const int& timestep);

std::string run_get_output(const std::string get_output, std::string data, int point);

int get_coeff_sc(const std::string& coeff_sc, 
	std::vector<double>& disp_x, 
	std::vector<double>& force0, 
	std::vector<double>& force1);

std::string run_gather_alya_output(const std::string get_alya_output, const int& run_id);

int parse_configfile(const std::string& config_file_name,
	std::string& nastin_dat,
	std::string& solidz_dat,
	std::string& create_data_point,
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
	double& rho_s_p2);

void get_output_data(const std::string get_output_sc, int& no_of_datapoints, std::vector<double>& disp_x, 
	std::vector<double>& force0, std::vector<double>& force1);

int save_coeff(const std::string file_name, const double& disp_x, const double& force0, const double& force1);

#endif /* HELPER_HPP_ */
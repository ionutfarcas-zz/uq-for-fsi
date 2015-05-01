#ifndef UQSIM_HPP_
#define UQSIM_HPP_

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <cstdlib>

class UQSimulation
{
public:
	virtual int data_decomp() const = 0;

	virtual std::vector<double> pre_processing() const = 0;
	virtual std::vector<double> pre_processing(const double& param1, const double& param2) const = 0;
	
	virtual void simulation(std::vector<double>& pre_proc_result) const = 0;

	virtual void post_processing() const = 0;

	virtual ~UQSimulation() {}
};

#endif /* UQSIM_HPP_ */
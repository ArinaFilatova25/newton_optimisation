#pragma once
#include "Function.h"
#include "Area.h"
#include <float.h>
#include <vector>
#include <random>
#include <chrono>
#include <sstream>
#include <cfloat>
#include "Stop_criterion.h"
#include "Optimization_method.h"
#include <Eigen/Dense>

class Newton_opt :
    public Optimization_method
{
private:
    std::vector<double> p;
    void correct_step_vector();
    std::pair<std::vector<double>, double> backtracking(std::vector<double> x);

public:
    Newton_opt();
    Newton_opt(Function* function, std::vector<double> x_0, Area area, Stop_criterion* stop_criterion);

    void optimization() override;
};





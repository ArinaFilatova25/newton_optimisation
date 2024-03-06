class Optimization_method;

#pragma once
#include "Function.h"
#include "Area.h"
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <sstream>

class Stop_criterion
{
protected:
    double eps;
    int max_num_of_iterations;
public:
    Stop_criterion();
    virtual ~Stop_criterion();
    Stop_criterion(double eps_, int max_num_of_iterations_);
    double get_eps();
    int get_max_num_of_iterations();

    virtual bool termination(Optimization_method* optimization_method) = 0;
};


// |f(x_{n+j}) − f(x_{n})| < eps, j = min{m: f(x_{n+m}) < f(x_{n})}
class Criterion_f_difference_min : public Stop_criterion {
public:
    Criterion_f_difference_min();
    explicit Criterion_f_difference_min(double eps_, int max_num_of_iterations_ = 100);
    bool termination(Optimization_method* optimization_method) override;
};

// число итераций больше n
class Criterion_max_iter : public Stop_criterion {
public:
    Criterion_max_iter();
    explicit Criterion_max_iter(int max_num_of_iterations_);
    bool termination(Optimization_method* optimization_method) override;
};

// число итераций с последнего улучшения больше n
class Criterion_num_iter_last_approx : public Stop_criterion {
public:
    Criterion_num_iter_last_approx();
    explicit Criterion_num_iter_last_approx(int max_num_of_iterations_);
    bool termination(Optimization_method* optimization_method) override;
};

 //  ||grad f(x_{n})|| < eps
class Criterion_grad_f : public Stop_criterion {
public:
    Criterion_grad_f();
    explicit Criterion_grad_f(double eps_, int max_num_of_iterations_ = 100);
    bool termination(Optimization_method* optimization_method) override;
};

// ||x_{n} - x_{n-1}|| < eps
class Criterion_x_difference : public Stop_criterion {
public:
    Criterion_x_difference();
    explicit Criterion_x_difference(double eps_, int max_num_of_iterations_ = 100);
    bool termination(Optimization_method* optimization_method) override;
};

// ||(f(x_{n}) - f(x_{n-1})/f(x_n)|| < eps
class Criterion_f_difference : public Stop_criterion {
public:
    Criterion_f_difference();
    explicit Criterion_f_difference(double eps_, int max_num_of_iterations = 100);
    bool termination(Optimization_method* optimization_method) override;
};


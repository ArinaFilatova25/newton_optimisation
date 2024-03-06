#include "Newton_opt.h"

Newton_opt::Newton_opt(){}

Newton_opt::Newton_opt(Function* function, std::vector<double> x_0, Area area,
    Stop_criterion* stop_criterion) : Optimization_method(function, x_0, area, stop_criterion) {
    p.resize(function->get_dim(), 0);
}

void Newton_opt::correct_step_vector() {
    std::vector<std::pair<double, double>> box = area.get_box(); //минимум и максимум для каждой координаты
    std::vector<double> x_n = seq_x_i.back();
    int dim = function->get_dim();

    // Проверка лежит ли вектор внутри области 
    bool vector_in_area = true;
    for (int i = 0; i < dim; ++i) {
        if (x_n[i] + p[i] < box[i].first || x_n[i] + p[i] > box[i].second) {
            vector_in_area = false;
        }
    }

    if (!vector_in_area) {
        for (int i = 0; i < dim; ++i) {
            p[i] = x_n[i] + p[i] < box[i].first ? box[i].first - x_n[i] : p[i];
            p[i] = x_n[i] + p[i] > box[i].second ? box[i].second - x_n[i] : p[i];
        }
    }
    else {
        double squar_sum = 0;
        for (int i = 0; i < dim; ++i) {
            squar_sum += p[i] * p[i];
        }

        double p_norm = std::sqrt(squar_sum);
        for (int i = 0; i < dim; ++i) {
            p[i] /= p_norm;
        }

        double min_alpha = DBL_MAX;
        double current_alpha = 0.0;
        for (int i = 0; i < dim; ++i) {
            current_alpha = (p[i] > 0) ? std::abs((box[i].second - x_n[i]) / p[i]) : std::abs((x_n[i] - box[i].first) / p[i]);
            min_alpha = (current_alpha < min_alpha) ? current_alpha : min_alpha;
        }

        for (int i = 0; i < dim; ++i) {
            p[i] *= min_alpha;
        }
    }
}


std::pair<std::vector<double>, double> Newton_opt::backtracking(std::vector<double> x) {
    double alpha = 1;
    double f_x = function->calculate(x);
    std::vector<double> new_x = x;
    int dim = function->get_dim();

    while (true) {
        for (int i = 0; i < dim; ++i) {
            new_x[i] = x[i] + alpha * p[i];
        }
        double new_f_x = function->calculate(new_x);

        if (new_f_x <= f_x) {
            return std::make_pair(new_x, new_f_x);
        }

        alpha /= 2.0; 
    }
}


void Newton_opt::optimization() {
   // num_of_iter = 0;
   // num_of_iter_since_last_approx = 0;
    int dim = function->get_dim();
    double eps = stop_criterion->get_eps();

    while (!stop_criterion->termination(this)) {
        std::vector<std::vector<double>> hess = function->hessian(seq_x_i.back(), eps / 10, area);

        Eigen::MatrixXd hessian_matrix(dim, dim);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                hessian_matrix(i, j) = hess[i][j];
            }
        }

        Eigen::MatrixXd inverse_hessian_matrix = hessian_matrix.inverse();

        std::vector<std::vector<double>> inverse_hessian(dim, std::vector<double>(dim, 0));
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                inverse_hessian[i][j] = inverse_hessian_matrix(i, j);
            }
        }

        std::vector<double> grad = function->gradient(seq_x_i.back(), eps / 10, area);

        std::vector<double> p_new(dim, 0);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                p_new[i] = -1 * inverse_hessian[i][j] * grad[j];
            }
        }
        p = p_new; 

        correct_step_vector();

        std::pair<std::vector<double>, double> new_x_f = backtracking(seq_x_i.back());
        seq_x_i.push_back(new_x_f.first);
        seq_f_i.push_back(new_x_f.second);

        ++num_of_iter;
        ++num_of_iter_since_last_approx;
    }
}

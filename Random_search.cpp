#include "Random_search.h"

Random_search::Random_search(Function* function, std::vector<double> x_0, Area area,
    Stop_criterion* stop_criterion, double p_, double delta_, double alpha_)
    : Optimization_method(function, x_0, area, stop_criterion), p(p_), delta(delta_), curr_delta(delta_), alpha(alpha_),
    seed(std::chrono::system_clock::now().time_since_epoch().count()), generator(seed), distribution(0, 1) {
    std::vector<std::pair<double, double>> box = area.get_box();
    int dimensions = function->get_dim();

    if (delta_ == -1) {
        delta = box[0].second - box[0].first;
        for (int i = 1; i < dimensions; ++i) {
            delta = box[i].second - box[i].first < delta ? box[i].second - box[i].first : delta;
        }
        delta /= 10;
        curr_delta = delta;
    }
}

void Random_search::optimization() {
    num_of_iter = 0;
    num_of_iter_since_last_approx = 0;

    int dimensions = function->get_dim();
    std::vector<std::pair<double, double>> box = area.get_box();
    bool is_in_small_area = false;

    while (!stop_criterion->termination(this)) {
        ++num_of_iter;
        ++num_of_iter_since_last_approx;
        std::vector<double> new_x;

        if (1 - p > distribution(generator)) { // Генерациия равномерного распределения во всем D.
            for (int i = 0; i < dimensions; ++i) {
                double min = box[i].first, max = box[i].second;
                new_x.push_back(distribution(generator) * (max - min) + min);
            }
            is_in_small_area = false;
        }
        else { // Генерация в пересечении D и B(x_n, delta)
            for (int i = 0; i < dimensions; ++i) {
                double min = seq_x_i.back()[i] - curr_delta > box[i].first ? seq_x_i.back()[i] - curr_delta : box[i].first,
                    max = seq_x_i.back()[i] + curr_delta < box[i].second ? seq_x_i.back()[i] + curr_delta : box[i].second;
                new_x.push_back(distribution(generator) * (max - min) + min);
            }
            is_in_small_area = true;
        }

        double new_f;
        if ((new_f = function->calculate(new_x)) < seq_f_i.back()) {

            if (is_in_small_area) {
                curr_delta *= alpha;
            }
            else {
                curr_delta = delta;
            }

            seq_f_i.push_back(new_f);
            seq_x_i.push_back(new_x);
            num_of_iter_since_last_approx = 0;
        }
    }
}

#include <iostream>
#include <memory>
#include <locale>
#include "Area.h"
#include "Function.h"
#include "Optimization_method.h"
#include "Stop_criterion.h"
#include "Newton_opt.h"
#include "Random_search.h"

enum FunctionType {
    FUNC_BUTA = 1,
    FUNC_ROSENBROCK = 2,
    FUNC_BILL = 3
};

Function* createFunction(int functionType) {
    switch (functionType) {
    case FUNC_BUTA:
        return new Function1();
    case FUNC_ROSENBROCK:
        return new Function2();
    case FUNC_BILL:
        return new Function3();
    default:
        return nullptr;
    }
}

int main() {
    setlocale(LC_ALL, "Rus");

    std::unique_ptr<Function> function;
    int function_nt = -1;

    std::cout << " Выберите функцию:\n"
        << " 1) Функция Бута (dim = 2);\n"
        << " 2) Функция Розенброка (dim = 3);\n"
        << " 3) Функция Билла (dim = 2).\n";

    std::cin >> function_nt;

    function = std::unique_ptr<Function>(createFunction(function_nt));

    int dim = function->get_dim();
    std::vector<double> x_0(dim);
    std::vector<std::pair<double, double>> box(dim);

    double coord = 0;
    try {
        std::cout << " Введите начальную точку:" << std::endl;
        for (int i = 0; i < dim; ++i) {
            std::cin >> coord;
            x_0[i] = coord;
        }

        std::cout << " Введите координаты границ области минимизации:" << std::endl;
        for (int i = 0; i < dim; ++i) {
            std::cout << " Ось " << i + 1 << ":" << std::endl;
            std::cin >> box[i].first >> box[i].second;
            if (box[i].first > x_0[i] || box[i].second < x_0[i]) {
                std::cout << " Начальная точка не находится в области минимизации." << std::endl;
                throw(-1);
            }
        }
    }
    catch (...) {
        return -1;
    }

    Area area(box);

    Optimization_method* optimization_method;
    Stop_criterion* stop_criterion;

        int optimization_method_int = -1, stop_criterion_int = -1, max_num_of_iter= -1;
        double eps = -1;
        std::cout << " Выберите метод оптимизации:\n"
            << " 1) Метод Ньютона (backtraking);\n"
            << " 2) Случайный поиск.\n";

        std::cin >> optimization_method_int;
    
        switch (optimization_method_int) {
        case 1:
            std::cout << " Выберите критерий остановки:\n"
                << " 1) ||grad f(x_{n})|| < eps;\n"
                << " 2) ||x_{n} - x_{n-1}|| < eps;\n"
                << " 3) ||(f(x_{n}) - f(x_{n-1}))/f(x_n)|| < eps.\n";

            std::cin >> stop_criterion_int;

            switch (stop_criterion_int) {
            case 1:
                std::cout << "Введите значение eps: ";
                std::cin >> eps;
                stop_criterion = new Criterion_grad_f(eps);
                break;
            case 2:
                std::cout << "Введите значение eps: ";
                std::cin >> eps;
                stop_criterion = new Criterion_x_difference(eps);
                break;
            case 3:
                std::cout << "Введите значение eps: ";
                std::cin >> eps;
                stop_criterion = new Criterion_f_difference(eps);
                break;
            default:
                stop_criterion = nullptr;
                break;
            }

            optimization_method = new Newton_opt(function.get(), x_0, area, stop_criterion);
            break;

        case 2:
            double p, delta , alpha ;
            std::cout << "Введите параметры p, delta и alpha: ";
            std::cin >> p >> delta >> alpha;
            std::cout << " Выберите критерий остановки для случайного поиска:\n"
                << " 1) |f(x_{n+j}) − f(x_{n})| < eps, j = min{m: f(x_{n+m}) < f(x_{n})};\n"
                << " 2) Число итераций больше n;\n"
                << " 3) Число итераций с последнего улучшения больше n.\n";

            std::cin >> optimization_method_int;

            switch (optimization_method_int) {
            case 1:
                double eps_f_difference_min;
                std::cout << "Введите значение eps: ";
                std::cin >> eps_f_difference_min;
                stop_criterion = new Criterion_f_difference_min(eps_f_difference_min);
                break;
            case 2:
                int max_num_of_iter;
                std::cout << "Введите значение n: ";
                std::cin >> max_num_of_iter;
                stop_criterion = new Criterion_max_iter(max_num_of_iter);
                break;
            case 3:
                int max_num_of_iter_last_approx;
                std::cout << "Введите значение n: ";
                std::cin >> max_num_of_iter_last_approx;
                stop_criterion = new Criterion_num_iter_last_approx(max_num_of_iter_last_approx);
                break;
            default:
                p = 0;
                delta = 0; 
                alpha = 0; 
                stop_criterion = nullptr;
                break;
            }

            optimization_method = new Random_search(function.get(), x_0, area, stop_criterion, p, delta, alpha);
            break;


        default:
            stop_criterion = new Criterion_grad_f(eps);
            optimization_method = new Newton_opt(function.get(), x_0, area, stop_criterion);
            break;
        }

    optimization_method->optimization();

    std::cout << "\nНаименьшее значение функции: " << optimization_method->get_seq_f_i().back() <<
        "\nТочка минимума: ";
    std::vector<double> x_n = optimization_method->get_seq_x_i().back();
    for (int i = 0; i < dim; ++i) {
        std::cout << x_n[i] << " ";
    }
    std::cout << "\nКоличество итераций: " << optimization_method->get_num_of_iter();

    return 0;
}



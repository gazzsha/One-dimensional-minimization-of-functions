#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#define esp 1.0e-4
#define t_ 4.0

double y(const double& x) { return 10 * x * log(x) - pow(x, 2) / 2; }
double fun(const double& x) { return x * atan(x) - 1.0 / 2.0 * log(1 + pow(x, 2)); }
double fun2(const double& x) { return exp(x) - 1 - x - pow(x, 2) / 2 - pow(x, 3) / 6; }
double derivative(double& x) {
    const double delta = 1.0e-3;
    double x1 = x - delta;
    double x2 = x + delta;
    double y1 = y(x1);
    double y2 = y(x2);
    return (y2 - y1) / (x2 - x1);
}

double derivate_second(double& x) {
    const double delta = 1.0e-3;
    double x1 = x - delta;
    double x2 = x + delta;
    double y1 = derivative(x1);
    double y2 = derivative(x2);
    return (y2 - y1) / (x2 - x1);
}

double df1(double& x) {
    const double delta = 1.0e-3;
    double x1 = x - delta;
    double x2 = x + delta;
    double y1 = fun(x1);
    double y2 = fun(x2);
    return (y2 - y1) / (x2 - x1);
}

double df2(double& x) {
    const double delta = 1.0e-3;
    double x1 = x - delta;
    double x2 = x + delta;
    double y1 = df1(x1);
    double y2 = df1(x2);
    return (y2 - y1) / (x2 - x1);
}
double Bitwise_search_method(double x_l, double& x_r, const double& t, int& count) {
    double delta = (x_r - x_l) / t;
    std::cout << std::fixed << std::setprecision(3) << "FIRST: " << std::setw(10) << "DELTA = " << delta
              << " "
              << "ESP = " << esp << '\n';
    std::cout << std::setw(4) << "X" << std::setw(16) << "Y(X)" << std::setw(16) << "Count" << '\n';
    do {
        count++;
        std::cout << std::fixed << std::setprecision(6) << x_l << std::setw(15) << y(x_l) << std::setw(15)
                  << count << '\n';
        x_l += delta;
    } while (y(x_l) <= y(x_l - delta) && x_l != x_r);
    count++;
    std::cout << std::fixed << std::setprecision(6) << x_l << std::setw(15) << y(x_l) << std::setw(15)
              << count << '\n';
    delta = delta / t;
    std::cout << std::fixed << std::setprecision(3) << "SECOND: " << std::setw(9) << "DELTA = " << delta
              << " "
              << "ESP = " << esp << '\n';

    std::cout << std::setw(4) << "X" << std::setw(16) << "Y(X)" << std::setw(16) << "Count" << '\n';
    do {
        count++;
        std::cout << std::setprecision(6) << x_l << std::setw(15) << y(x_l) << std::setw(15) << count << '\n';
        x_l -= delta;
    } while (y(x_l + delta) >= y(x_l) && x_l != x_r);
    count++;
    std::cout << std::fixed << std::setprecision(6) << x_l << std::setw(15) << y(x_l) << std::setw(15)
              << count << '\n';
    x_l += delta;
    return x_l;
}

void print_midpoint(const double& x_l, const double& x_r, const double& x_m, const double& der,
                    const int& count) {
    std::cout << std::fixed << std::setprecision(6) << x_l << std::setw(9) << x_r << std::setw(12) << x_m
              << std::setw(12) << der << std::setw(17) << count << '\n';
}

double midpoint_method(const double& x_l, const double& x_r, int& count) {
    count++;
    double x_mid = (x_r + x_l) / 2;
    double der = (derivative(x_mid));
    print_midpoint(x_l, x_r, x_mid, der, count);
    while (fabs(der) > esp) {
        if (der > 0)
            return midpoint_method(x_l, x_mid, count);
        else
            return midpoint_method(x_mid, x_r, count);
    }
    return x_mid;
}

void print_Newton_Raphson_method(const double& x, const double& x_next, const double& tay, const double& zap,
                                 const int& count) {
    std::cout << std::fixed << std::setprecision(4) << x << std::setw(9) << x_next << std::setw(12) << tay
              << std::setw(12) << zap << std::setw(12) << count << '\n';
}
double Newton_Raphson_method(double& x_l, int& count) {
    count++;
    double x_del = x_l - derivative(x_l) / derivate_second(x_l);
    double t_opt = pow(derivative(x_l), 2) / (pow(derivative(x_l), 2) + pow(derivative(x_del), 2));
    double x_next = x_l - t_opt * derivative(x_l) / derivate_second(x_l);
    print_Newton_Raphson_method(x_l, x_next, t_opt, x_next - x_l, count);
    if (x_next - x_l > esp) {
        return Newton_Raphson_method(x_next, count);
    }
    return x_next;
}

double Newton_Raphson_method_(double& x_l, int& count) {
    count++;
    double x_del = x_l - df1(x_l) / df2(x_l);
    double t_opt = pow(df1(x_l), 2) / (pow(df1(x_l), 2) + pow(df1(x_del), 2));
    double x_next = x_l - t_opt * df1(x_l) / df2(x_l);
    print_Newton_Raphson_method(x_l, x_next, t_opt, x_next - x_l, count);
    if (x_next - x_l > esp) {
        return Newton_Raphson_method_(x_next, count);
    }
    return x_next;
}

void print_fibonacci_method(const double& xl, const double& xr, const double& l, const double& m,
                            const double& y1, const double& y2, const int& count) {
    std::cout << std::fixed << std::setprecision(6) << xl << std::setw(9) << xr << std::setw(9) << l
              << std::setw(9) << m << std::setw(12) << y1 << std::setw(12) << y2 << std::setw(12) << count
              << '\n';
}
double Fibo(int n) {
    if (n == 1 || n == 2 || n == 0) return 1;
    return Fibo(n - 1) + Fibo(n - 2);
}
double fibonacci_method(double x_l, double x_r, int& count) {
    count = 0;
    double x_ = 0;
    int n = ((x_r - x_l) / esp);
    int k = 1;
    while (Fibo(k) < n) {
        k++;
    }
    n = k;
    double x1, x2, y1, y2;
    for (int i = 0; i < n - 1; ++i) {
        x1 = x_l + Fibo(n - i - 1) / Fibo(n + 1 - i) * (x_r - x_l);
        x2 = x_l + Fibo(n - i) / Fibo(n - i + 1) * (x_r - x_l);
        y1 = y(x1);
        y2 = y(x2);
        count++;
        print_fibonacci_method(x_l, x_r, x1, x2, y1, y2, count);
        if (y1 > y2) {
            x_l = x1;
            ;
        } else if (y2 > y1) {
            x_r = x2;
        }
    }
    return (x_l + x_r) / 2;
}
void print_parabola_method(const double& x, const double& y, const double del, const int& count) {
    std::cout << std::fixed << std::setprecision(6) << x << std::setw(15) << y << std::setw(15) << del
              << std::setw(9) << count << '\n';
}
double parabola_method(double x1, double x2, double x3, int& count) {
    double y1, y2, y3, x_, y_, a0, a1, a2;
    double e = 1.0e-6;
    y1 = fun2(x1);
    y2 = fun2(x2);
    y3 = fun2(x3);
    a0 = y1;
    a1 = (y2 - y1) / (x2 - x1);
    a2 = 1 / (x3 - x2) * ((y3 - y1) / (x3 - x1) - (y2 - y1) / (x2 - x1));
    x_ = 0.5 * (x1 + x2 - a1 / a2);
    y_ = fun2(x_);
    count++;
    print_parabola_method(x_, y_, x_ - x2, count);
    while (fabs(x_ - x2) > e) {
        if (y_ <= y2) {
            if (x_ < x2) {
                x3 = x2;
                x2 = x_;
            } else {
                x1 = x2;
                x2 = x_;
            }
        } else {
            if (x_ < x2) {
                x1 = x_;
            } else {
                x1 = x2;
                x2 = x_;
            }
        }
        y1 = fun2(x1);
        y2 = fun2(x2);
        y3 = fun2(x3);
        a1 = (y2 - y1) / (x2 - x1);
        a2 = 1 / (x3 - x1) * ((y3 - y1) / (x3 - x1) - (y2 - y1) / (x2 - x1));
        x_ = 0.5 * (x1 + x2 - a1 / a2);
        y_ = fun2(x_);
        count++;
        print_parabola_method(x_, y_, x_ - x2, count);
    }

    return x_;
}

void initial_approximation(double& a, double& b) {
    double delta = 0.1;
    while ((df1(a) >= 0 && df2(a) >= 0) || (df1(a) <= 0 && df2(a) <= 0)) {
        a -= delta;
    }
    while ((df1(b) >= 0 && df2(b) >= 0) || (df1(b) <= 0 && df2(b) <= 0)) {
        b += delta;
    }
}

void print_find_method(const double& x, const double& y, const int& n) {
    std::cout << std::fixed << std::setprecision(6) << x << std::setw(15) << y << std::setw(17) << n << '\n';
}
double find_method(const double& a, const double& b, int& count) {
    count = 0;
    double x0 = 0;
    double x1 = 0;
    double x;
    double y1, y2;
    double sh = (b - a) / 4;
    x0 = a;
    y1 = y(x0);
    count++;
label:
    print_find_method(x0, y1, count);
    count++;
    x1 = x0 + sh;
    y2 = y(x1);
    if (y1 >= y2) {
        x0 = x1;
        y1 = y2;
        if (x0 >= a && x0 <= b) goto label;
    }
    if (fabs(sh) < esp) {
        x = x0;
    } else {
        sh = -sh / 4;
        x0 = x1;
        y1 = y2;
        goto label;
    }
    return x;
}

int main() {
    double xl = 0.1;
    double xr = 1;
    int count = 0;
    // firsr
    std::cout << "First Method:" << '\n';
    std::cout << std::fixed << std::setprecision(3) << "X" << std::setw(15) << "Y" << std::setw(24) << "Count"
              << '\n';
    double x_min_find = find_method(xl, xr, count);
    std::cout << "X(min) = " << x_min_find << " "
              << "Y(min) = " << y(x_min_find) << '\n';

    // double x_min = Bitwise_search_method(xl, xr, t_, count);
    // double y_min = y(x_min);
    // std ::cout << "RESULT: ";
    // std::cout << "Fisrt Method:" << '\n';
    // std::cout << "X(min) = " << x_min << " "
    //           << "Y(min) = " << y_min << '\n';

    // фибо
    std::cout << "Second Method:" << '\n';
    std::cout << std::fixed << std::setprecision(3) << "X(Left)" << std::setw(10) << "X(Right)"
              << std::setw(2) << "l" << std::setw(9) << "M" << std::setw(15) << "F(l)" << std::setw(12)
              << "F(M)" << std::setw(14) << "Count" << '\n';
    double x_min_Fibo = fibonacci_method(xl, xr, count);
    std::cout << "X(min) = " << x_min_Fibo << " "
              << "Y(min) = " << y(x_min_Fibo) << '\n';
    // 3
    std::cout << "Third Method:" << '\n';
    std::cout << std::fixed << std::setprecision(3) << "X(Left)" << std::setw(10) << "X(Right)"
              << std::setw(10) << "X(mid)" << std::setw(15) << "Derivative" << std::setw(14) << "Count"
              << '\n';
    count = 0;
    double x_min_midpoint = midpoint_method(0.1, 1.0, count);
    std::cout << "X(min) = " << x_min_midpoint << " "
              << "Y(min) = " << y(x_min_midpoint) << '\n';

    // 4
    std::cout << "Fourth Method:" << '\n';
    std::cout << std::fixed << std::setprecision(3) << "X(Cur.)" << std::setw(9) << "X(Next)" << std::setw(9)
              << "T(k)" << std::setw(13) << "Delta" << std::setw(13) << "Count" << '\n';
    count = 0;
    double x_min_Newton_Raphson_method = Newton_Raphson_method(xl, count);
    std::cout << "X(min) = " << x_min_Newton_Raphson_method << " "
              << "Y(min) = " << y(x_min_Newton_Raphson_method) << '\n';

    // parabola
    std::cout << "Parabola Method:" << '\n';
    printf("%s", "x1 = -5, x2=0.4, x3=5\n");
    std::cout << std::fixed << std::setprecision(3) << "X" << std::setw(15) << "Y" << std::setw(19) << "Delta"
              << std::setw(12) << "Count" << '\n';
    count = 0;
    double x_min_parabola = parabola_method(-5, 0.4, 5, count);
    std::cout << "X(min) = " << x_min_parabola << " "
              << "Y(min) = " << fun2(x_min_parabola) << '\n';
    // Newton
    printf("%s", "Началльное приближение:\n");
    double a = 0;
    double b = 0;
    initial_approximation(a, b);
    printf("%s %lf %s %lf %s", "Диапозон: [", a, ";", b, "]\n");
    std::cout << std::fixed << std::setprecision(3) << "X(Cur.)" << std::setw(9) << "X(Next)" << std::setw(9)
              << "T(k)" << std::setw(13) << "Delta" << std::setw(13) << "Count" << '\n';
    count = 0;
    double x_min_Newton_Raphson_method_ = Newton_Raphson_method_(a, count);
    std::cout << "X(min) = " << x_min_Newton_Raphson_method_ << " "
              << "Y(min) = " << fun(x_min_Newton_Raphson_method_) << '\n';

    return 0;
}
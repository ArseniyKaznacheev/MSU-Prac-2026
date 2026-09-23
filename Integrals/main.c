#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern double f1(double x);
extern double f2(double x);
extern double f3(double x);
extern double f1p(double x);
extern double f2p(double x);
extern double f3p(double x);

typedef double afunc(double);

static unsigned long root_iterations = 0;

static double root(afunc* f, afunc* g, afunc* df, afunc* dg,
                                            double a, double b, double eps1)
{
    double Fa = f(a) - g(a);
    double Fb = f(b) - g(b);
    double mid = (a + b) / 2;
    double Fmid = f(mid) - g(mid);


    /*
    Поскольку производные F'(x) и F"(x) на отрезке [a, b] не меняют знак, 
    для проверки первого условия достаточно сравнить F(a) с 0 
    (при F(a) < 0 функция возрастает)

    0 != (Fa < 0.0) => f ^ (F'(x) > 0)     (Fmid < (Fa + Fb) / 2)) == 0 => выпукла вверх(F''(x) < 0)
    0 == (Fa > 0.0) => f v (F'(x) < 0)     (Fmid < (Fa + Fb) / 2)) != 0 => выпукла вниз(F''(x) > 0)

    */
    /*
    В методе касательных проводится касательная к кривой y = F(x) 
    в точке(b, F(b))в случае 1   (F'(x)F"(x) > 0)
    в точке(a, F(a)) в случае 2
    */

    double x = ((Fa < 0.0) == (Fmid < (Fa + Fb) / 2)) ? b: a;
   
    double x_new = x;
    root_iterations = 0;
    do {
        x = x_new;
        double Fx = f(x) - g(x);
        double dFx = df(x) - dg(x);
        x_new = x - Fx / dFx;
        root_iterations++;
    } while (fabs(x_new - x) > eps1);

    return x_new;
}

static double rect_sum(afunc* f, double a, double b, int n)
{
    // интеграл по заданному разбиению
    double h = (b - a) / n;
    double s = 0.0;
    for (int i = 0; i < n; i++)
        s += f(a + (i + 0.5) * h);
    return s * h;
}

static double integral(afunc* f, double a, double b, double eps2)
{
    int n = 10;
    double I_prev = rect_sum(f, a, b, n);
    double I_curr = I_prev;
    for (int k = 0; k < 25; k++) {
        n *= 2;
        I_curr = rect_sum(f, a, b, n);
        if (fabs(I_curr - I_prev) / 3 < eps2) {
            return I_curr;
        }
        I_prev = I_curr;
    }
    return I_curr;

    /*
        берется некоторое начальное число разбиений n0 = 10  и последовательно вычисляются значения
        In при n, равном 2n0, 4n0, 8n0 и т.д. 
        Известно правило Рунге:
                | I – In | ≅ p | In – I2n |
         p = 1/3
         Согласно этому правилу, когда на очередном шаге величина p*|In − I2n| < ε
            в качестве приближенного значения для I можно взять I2n
    */
}

//функции для площади
static double h1(double x) { return f1(x) - f3(x); }
static double h2(double x) { return f1(x) - f2(x); }

//тестовые функции
static double tf1(double x) {return x * x;}
static double tf1p(double x) {return 2.0 * x;}

static double tf2(double x) {return 2.0 - x;}
static double tf2p(double x) { (void)x;  return -1.0; }

static double tf3(double x) { return sin(x);}
static double tf3p(double x) {return cos(x);}

static double tf4(double x) {return - x * x * x;}
static double tf4p(double x) {return - 2.0 * x * x; }

static double tf5(double x) {return exp(x);}
static double tf5p(double x) {return exp(x);}

static double tf6(double x) { (void)x; return exp(1);}
static double tf6p(double x) { (void)x;  return 0.0; }

static afunc* tf_mas[] = {tf1, tf2, tf3, tf4, tf5, tf6};
static afunc* tfp_mas[] = {tf1p, tf2p, tf3p, tf4p, tf5p, tf6p};

static void print_help(void)
{
    for (int i = 0; i < 40; i++) {
        printf("-");
    }
    printf("\nUsage: integral [FLAGS][OPTIONS]\n");
    printf("  -h, --help                      print this help\n");
    printf("  -f, --functions                 available testing functions\n");
    printf("  -r, --root                      print abscissas of intersections\n");
    printf("  -i, --iterations                print count of iterationsof  root search\n");
    printf("  -R, --test-root F1:F2:A:B:E:R   test root function\n");
    printf("  -I, --test-integral F:A:B:E:R   test integral function\n\n");
    printf("F      - function number;\n");
    printf("[A; B] - Search segment;\n");
    printf("E      - epsilon;\n");
    printf("R      - correct answer (calculated analytically).\n");
    for (int i = 0; i < 40; i++) {
        printf("-");
    }
    printf("\n");
}

static void print_test_func(void) {
    printf("\nAvailable functions:\n");
    printf("1: y = x^2\n");
    printf("2: y = 2-x\n");
    printf("3: y = sin(x)\n");
    printf("4: y = -x^3\n");
    printf("5: y = exp(x)\n");
    printf("6: y = e\n\n");
}

static int run_test_root(const char* arg)
{
    int f1, f2;
    double a, b, eps, expected;
    if (sscanf_s(arg, "%d:%d:%lf:%lf:%lf:%lf", &f1, &f2, &a, &b, &eps, &expected) != 6) {
        fprintf(stderr, "Bad '--test-root' argument '%s'\n", arg);
        return 1;
    }
    if (f1 < 1 || f1 > 6 || f2 < 1 || f2 > 6) {
        fprintf(stderr, "No such function\n");
        return 1;
    }
    f1--;
    f2--;
    double x = root(tf_mas[f1], tf_mas[f2], tfp_mas[f1], tfp_mas[f2], a, b, eps);
    double abs_err = fabs(x - expected);
    double rel_err = (expected != 0.0) ? (abs_err / fabs(expected)) : abs_err;

    printf("\n%.7f %.7f %.10f\n", x, abs_err, rel_err);
    return 0;
}

static int run_test_integral(const char* arg)
{
    int f;
    double a, b, eps, expected;
    if (sscanf_s(arg, "%d:%lf:%lf:%lf:%lf", &f, &a, &b, &eps, &expected) != 5) {
        fprintf(stderr, "Bad '--test-integral' argument '%s'\n", arg);
        return 1;
    }
    if (f < 1 || f > 6) {
        fprintf(stderr, "No such function\n");
        return 1;
    }
    f--;
    double I = integral(tf_mas[f], a, b, eps);
    double abs_err = fabs(I - expected);
    double rel_err = (expected != 0.0) ? (abs_err / fabs(expected)) : abs_err;

    printf("%.7f %.7f %.10f\n", I, abs_err, rel_err);
    return 0;
}

int main(int argc, char* argv[])
{
    //printf("%.6f, %.6f, %.6f\n", f1(-2.5222234238), f2(0.6505195252), f2(1.2793531911));

    int flag_help = 0, flag_root = 0, flag_iter = 0, flag_func = 0;
    const char* test_root_arg = NULL;
    const char* test_int_arg = NULL;

    //парсинг командной строки  
    for (int i = 1; i < argc; i++) {
        const char* a = argv[i];

        if (a[0] == '-' && a[1] == '-') {
            if (!strcmp(a, "--help")) {
                flag_help = 1; continue;
            }
            if (!strcmp(a, "--root")) {
                flag_root = 1; continue;
            }
            if (!strcmp(a, "--iterations")) {
                flag_iter = 1; continue;
            }
            if (!strcmp(a, "--functions")) {
                flag_func = 1; continue;
            }
            if (!strcmp(a, "--test-root")) {
                if (++i >= argc) {
                    fprintf(stderr, "Option '--test-root' requires an argument\n");
                    return 1;
                }
                test_root_arg = argv[i];
                continue;
            }
            if (!strcmp(a, "--test-integral")) {
                if (++i >= argc) {
                    fprintf(stderr, "Option '--test-integral' requires an argument\n");
                    return 1;
                }
                test_int_arg = argv[i];
                continue;
            }

            fprintf(stderr, "Unknown option '%s'\n", a);
            return 1;
        }

        // короткие опции
        if (a[0] == '-' && a[1] != '\0' && a[2] == '\0') {
            switch (a[1]) {
            case 'h': flag_help = 1; break;
            case 'r': flag_root = 1; break;
            case 'i': flag_iter = 1; break;
            case 'f': flag_func = 1; break;
            case 'R':
                if (++i >= argc) {
                    fprintf(stderr, "Option '-R' requires an argument\n");
                    return 1;
                }
                test_root_arg = argv[i];
                break;
            case 'I':
                if (++i >= argc) {
                    fprintf(stderr, "Option '-I' requires an argument\n");
                    return 1;
                }
                test_int_arg = argv[i];
                break;
            default:
                fprintf(stderr, "Unknown option '%c'\n", a[1]);
                return 1;
            }
            continue;
        }
        fprintf(stderr, "Unexpected argument '%s'\n", a);
        return 1;
    }

    if (flag_help) {
        print_help();
    }
    if (flag_func) {
        print_test_func();
    }
    if (test_root_arg) {
        if (run_test_root(test_root_arg)) {
            return 1;
        }
    }
    if (test_int_arg) {
        if (run_test_integral(test_int_arg)) {
            return 1;
        }
    }

    // подсчет площади
    const double eps1 = 1e-4;
    const double eps2 = 1e-4;

    double xA = root(f1, f3, f1p, f3p, -3.0, -2.0, eps1);
    unsigned long iter_A = root_iterations;

    double xB = root(f2, f3, f2p, f3p, 0.0, 1.0, eps1);
    unsigned long iter_B = root_iterations;

    double xC = root(f1, f2, f1p, f2p, 1.0, 2.0, eps1);
    unsigned long iter_C = root_iterations;

    if (flag_root) {
        printf("\nIntersection points:\n");
        printf("A(f1 n f3): x = %.10f; y = %.6f\n", xA, f1(xA));
        printf("B(f2 n f3): x = %.10f; y = %.6f\n", xB, f2(xB));
        printf("C(f1 n f2): x = %.10f; y = %.6f\n", xC, f2(xC));
    }
    if (flag_iter) {
        printf("\nCount of root search iterations:\n");
        printf("A(f1 n f3): %lu\n", iter_A);
        printf("B(f2 n f3): %lu\n", iter_B);
        printf("C(f1 n f2): %lu\n", iter_C);
    }

    double S1 = integral(h1, xA, xB, eps2); //f1(x) - f3(x)
    double S2 = integral(h2, xB, xC, eps2); //f1(x) - f2(x)
    double S = S1 + S2;

    printf("\nS = %.6f\n", S);
    return 0;
}

/*
 ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⡔⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠜⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡎⠀C⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠂⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡰⠁⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢰⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠞⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠂⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠜⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠔⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⡠⠃⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⢀⡠⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⢀⠔⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⡠⠂⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡴⢼⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡈⠀⠀f2⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀f1⠀⢀⡠⠔⠉⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡠⠔⠊⠁⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠒⠦⠤⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⠤⠔⠊⠁⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠈⠁⠒⠠⠄⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⠠⠄⠒⠋⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⠒⠢⠄⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⡀⠤⠔⠒⠊⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⣁⣒⣠⣄⣠⡤⠤⠀⠐⠒⠂⠈⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠤⠤⠤⠀⠀⠀⠒⠒⠒⠒⠒⠈⠉⠉⠉⠉⠉⠀⠀⠀⠀A⠉⠙⠒⠤⠄⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀ ⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
                            ⠉⠛⠩⢍⣉                          ⢸              ⢹
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠑⠒⠤⠄⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠁⠒⠤⢄⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀f3⠀⠀⠀⠉⠑⠒⠤⢄⣀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡘⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠑⠒⠤⢄⣀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠰⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⢹⠒⠤⢄⣀⠀⠀⠀⠀⠀⠀⠀⠠⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠉⠑⠒⠤⢄⣀⡐⠁⠀B⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣸⣀⣀⣀⣀⣀⣀⣀⣤⣐⣉⣉⣑⣒⣤⣄⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣲⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢐⡂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠀⠀⠀⢀⠔⠊⠁⠀⠀⠀⠀⠀⢖⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡆⠉⠑⠒⠤⢄⡀⠀⠀⠀⠀⠀⠀⠀⠀⣲⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠔⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠑⠒⠤⢄⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡌⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠑⠂⠤⢀⡀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡘⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠐⠒⠤⢄⣀
 */
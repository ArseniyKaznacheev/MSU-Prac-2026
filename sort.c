#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include <float.h>

unsigned int com = 0, trans = 0;

static void print_mas(double* mas, unsigned int length) {
    for (unsigned int i = 0; i < length; i++) {
        printf("%e ", mas[i]);
    }
    printf("\n");
}


static int compar(double x, double y) {
    com++;
    if (fabs(x) > fabs(y)) {
        return 1;
    }
    else {
        return 0;
    }
}
static int test(double* mas, double* mas_sort, unsigned int length) {
    // проверяем отсортированность

    for (unsigned int i = 0; i < (length - 1); i++) {
        if (compar(mas_sort[i + 1], mas_sort[i])) {
            return 1;
        }
    }
    com = 0;
    trans = 0;
    //проверяем биекцию

    unsigned int check_length = length;
    for (unsigned int i = 0; i < length; i++) {
        unsigned int j = 0;
        while ((j < check_length) &&(mas[i] != mas_sort[j])) {
            j++;
        }
        if (j == check_length) {
            return 1;
        }

        for (unsigned int k = j; k < (check_length - 1); k++) {
            mas_sort[k] = mas_sort[k + 1];
        }
        mas_sort[check_length - 1] = 0;
        check_length--;
    }

    return 0;
}

static void sort(int L, int R, double* mas) {
    int i = L;
    int j = R;
    double a, b;

    a = mas[(L + R) / 2];
    do {
        while (compar(mas[i], a)) {
            i++;
        }
        while (compar(a, mas[j])) {
            j--;
        }

        if (i <= j) {
            b = mas[i];
            mas[i] = mas[j];
            mas[j] = b;
            trans++;
            i++;
            j--;
        }
    } while (i <= j);

    if (L < j) {
        sort(L, j, mas);
    }
    if (i < R) {
        sort(i, R, mas);
    }
}
static void q_sort_rec(double* mas, unsigned int length) {
    if (length > 1) {
        sort(0, length - 1, mas);
    }
}

static void bubblesort(double* mas, unsigned int length) {
    for (unsigned int i = 0; i < length; i++) {
        for (unsigned int j = length - 1; j > i; j--) {
            if (compar(mas[j], mas[j-1])) {
                double temp = mas[j];
                mas[j] = mas[j - 1];
                mas[j - 1] = temp;
                trans++;
            }
        }
    }
}

/*
static void shell_sort(double* mas, unsigned int length) {
    for (unsigned int s = (length / 2); s > 0; s /= 2) { //шаг
        for (unsigned int i = s; i < length; ++i) {
            for (unsigned int j = i; j >= s && compar(mas[j], mas[j - s]); j -= s) { // сортировка вставками
                double temp = mas[j];
                mas[j] = mas[j - s];
                mas[j - s] = temp;
                trans++;
            }
        }
    }
}
*/
static int generate_rand(double* mas, unsigned int length) {
    if (mas == NULL || length == 0) {
        return 1;
    }
    //srand((unsigned int)time(NULL));
    for (unsigned int i = 0; i < length; i++) {
        double f = (double)rand() / RAND_MAX;
        f = (2.0 * f - 1.0) * DBL_MAX;
        //f = DBL_MIN + f * DBL_MAX - f * DBL_MIN;
        mas[i] = f;
    }
    return 0;
}

static int generate_sorted(double* mas, unsigned int length) {
    if (mas == NULL || length == 0) {
        return 1;
    }
    //srand((unsigned int)time(NULL));
    int com = 0, trans = 0;
    for (unsigned int i = 0; i < length; i++) {
        double num = (double)rand() / RAND_MAX;
        num = (2.0 * num - 1.0) * DBL_MAX;
        //num = (rand() % 2 == 0) ? num : -num;

        int j = i;
        while (j > 0) {
            if (compar(num, mas[j - 1])) {
                mas[j] = mas[j - 1];
                j--;
            }
            else {
                break;
            }
        }
        mas[j] = num;
    }
    return 0;
}

static int generate_sorted_rev(double* mas, unsigned int length) {
    if (mas == NULL || length == 0) {
        return 1;
    }
    //srand((unsigned int)time(NULL));
    int com = 0, trans = 0;
    for (unsigned int i = 0; i < length; i++) {
        double num = (double)rand() / RAND_MAX;
        num = (2.0 * num - 1.0) * DBL_MAX;
        //num = (rand() % 2 == 0) ? num : -num;

        int j = i;
        while (j > 0) {
            if (compar(num, mas[j - 1])) {
                mas[j] = mas[j - 1];
                j--;
            }
            else {
                break;
            }
        }
        mas[j] = num;
    }
    for (unsigned int i = 0; i < (length / 2); i++) {
        double temp = mas[i];
        mas[i] = mas[length - 1 - i];
        mas[length - 1 - i] = temp;
    }
    return 0;
}

static void print_res(unsigned int length, int* mas, char* sort_type, char* param) {

    printf("                    %s, %s\n", sort_type, param);
    for (int i = 0; i < 62; i++) {
        printf("-");
    }
    printf("\n");

    printf("|                         Mass Number               Mean     |\n|");
    for (int i = 0; i < 10; i++) {
        printf(" ");
    }
    for (int i = 0; i < 40; i++) {
        printf("-");
    }
    for (int i = 0; i < 10; i++) {
        printf(" ");
    }
    printf("|\n|");

    printf("n         1         2         3         4                   |\n|");
    printf("%-10u", length);
    int sum = 0;
    for (int j = 0; j < 4; j++) {
        printf("%-10d", mas[j]);
        sum += mas[j];
    }
    //printf("   ");
    printf("%-10.2f|\n", (double)sum / 4);

    for (int i = 0; i < 62; i++) {
        printf("-");
    }
    printf("\n\n\n");
}


int main(void) {
    
    srand((unsigned int)time(NULL));
    /*
    double f = (double)rand() / RAND_MAX;
    f = (2.0 * f - 1.0) * DBL_MAX;
    printf("%e\n", f);
    */
    unsigned int length;
    unsigned int com_mas_1[4]; // массив количеств сравнений для Qsort
    unsigned int trans_mas_1[4];
    unsigned int com_mas_2[4]; // для Пузырька
    unsigned int trans_mas_2[4];

    unsigned int logs_qsort[4];
    unsigned int logs_bubble[4];

    if (scanf("%u", &length) != 1) {
        printf("Reading array length error");
        return 1;
    }
    if (!length) {
        printf("Zero length");
        return 0;
    }

    double* mas = malloc(sizeof(double) * length);
    double* mas_copy = malloc(sizeof(double) * length);

    if (mas == NULL || mas_copy == NULL) {
        printf("Mem error");
        return 0;
    }

    int (*gen_func[])(double*, unsigned int) = {generate_sorted, generate_sorted_rev, generate_rand, generate_rand};


    for (int i = 0; i < 4; i++) {
        if (gen_func[i](mas, length)) {
            printf("Generation error\n");
            return 1;
        }
        com = 0;
        trans = 0;

        memcpy(mas_copy, mas, sizeof(double) * length);

        q_sort_rec(mas_copy, length);
        //print_mas(mas_copy, length);//
        //print_mas(mas, length);//

        com_mas_1[i] = com;
        trans_mas_1[i] = trans;
        com = 0;
        trans = 0;
        logs_qsort[i] = test(mas, mas_copy, length);

        memcpy(mas_copy, mas, sizeof(double) * length);
        //shell_sort(mas_copy, length);
        bubblesort(mas_copy, length);
        //print_mas(mas_copy, length);//
        //print_mas(mas, length);

        com_mas_2[i] = com;
        trans_mas_2[i] = trans;
        com = 0;
        trans = 0;
        logs_bubble[i] = test(mas, mas_copy, length);
        //printf("\n");
        //printf("\n");
    }

    free(mas);
    free(mas_copy);
    char* param1 = "comparisons count";
    char* param2 = "exchanges count";
    char* sort_type2 = "BubbleSort";
    char* sort_type1 = "Qsort";

    print_res(length, com_mas_1, sort_type1, param1);
    print_res(length, trans_mas_1, sort_type1, param2);

    print_res(length, com_mas_2, sort_type2, param1);
    print_res(length, trans_mas_2, sort_type2, param2);

    printf("tests_Qsort:\n");
    for (int i = 0; i < 4; i++) {
        printf("%d: %s", i+1, !logs_qsort[i] ? "OK": "Not OK");
        printf("\n");
    }

    printf("tests_BubbleSort: \n");
    for (int i = 0; i < 4; i++) {
        printf("%d: %s ", i+1, !logs_bubble[i] ? "OK" : "Not OK");
        printf("\n");
    }

    return 0;
}
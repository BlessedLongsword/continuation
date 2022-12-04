/*Antonio Tigri - Mètodes Numèrics (Pràctica 2: Continuació)*/

/* Includes */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* Defines */
#define NUM 30000
#define DELTA 0.005
#define TOL 1e-5
#define PRE 1e-8
#define ITER 5

/* Header Functions */
void generate_curve(double x0, double y0, const int niub[], int num, double delta, double tol, double pre, int iter,
                    int prints, FILE * fptr);
double f(double x, double y, const int niub[]);
double dx_f(double x, double y, const int niub[]);
double dy_f(double x, double y, const int niub[]);
double h(double x, double y, double x0, double y0, const int niub[]);
double dx_h(double x, double y, const int niub[]);
double dy_h(double x, double y, const int niub[]);
double dist(double x, double y, double p, double q);
double dx_dist(double x, double y, double p, double q);
double dy_dist(double x, double y, double p, double q);
double compute_b(const int niub[], int comp);
double compute_initial_points(const int niub[], const double distances[], double const b[],
                              int comp, int point, int dist_idx);

/* ------------------------- MAIN ------------------------- */

int main() {
    int i, j, counter = 0;
    int prints;
    int  niub[8];
    double init_points[16][2];
    double distances[4] = {0.2, 0.1, -0.05, -0.15};
    FILE * fptr;

    printf("\nInput your (niub) number with spaces (i.e. 2 0 0 3 4 5 6 1):\n");
    for (i = 0; i < 8; i++) {
        scanf("%d", &niub[i]);
    }

    printf("\nThe Points Pi are: ");
    for (i = 0; i < 4; i++) {
        printf("P%d = (%d, %d) ", i, niub[2*i], niub[2*i+1]);
    }

    printf("\n\nComputing B...\n");
    double b[2] = {compute_b(niub, 0), compute_b(niub, 1)};
    printf("(B = %1.3f, %1.3f)\n", b[0], b[1]);

    printf("\nComputing Initial Points...\n");
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            init_points[counter][0] = compute_initial_points(niub, distances, b, 0, i, j);
            init_points[counter][1] = compute_initial_points(niub, distances, b, 1, i, j);
            counter += 1;
        }
    }
    printf("Initial Points: [");
    for(i = 0; i < 15; i++) {
        printf("(%1.3f, %1.3f), ", init_points[i][0], init_points[i][1]);
    } printf("(%1.3f, %1.3f)]\n\n", init_points[15][0], init_points[15][1]);

    fptr = fopen("curves.txt", "w");

    do {
        printf("Do you want to print predictions and corrections? [0 = no, 1 = yes]\n");
        scanf("%d", &prints);
        if (prints < 0 || prints > 1) {
            printf("Invalid input! ");
        }
    } while (prints < 0 || prints > 1);

    for(i = 0; i < 16; i++) {
        printf("\nGenerating curve #%d:\n", i+1);
        generate_curve(init_points[i][0], init_points[i][1], niub, NUM, DELTA, TOL, PRE, ITER, prints, fptr);
        printf("#Starting point reached, finished curve #%d.\n\n", i+1);
        fprintf(fptr, "\n#Starting point reached, finished curve #%d.\n\n", i+1);
    }

    fclose(fptr);

    return 0;
}

/* ----------------- Main Functions -----------------*/
void generate_curve(double x0, double y0, const int niub[], int num, double delta, double tol, double pre, int iter,
                    int prints, FILE * fptr) {
    int iter_counter = 0;
    double x_init = x0, y_init = y0;
    double x_pred, y_pred;
    double tang_x, tang_y;

    do {

        /* Prediction */
        tang_x = dy_h(x_init, y_init, niub);
        tang_y = -dx_h(x_init, y_init, niub);

        if (dist(tang_x, tang_y, 0, 0) < tol) {
            printf("\nThe gradient is too small! Stopping curve generation at iteration %d\n\n", iter_counter);
            break;
        } else {
            x_pred = x_init + delta * (tang_x / dist(tang_x, tang_y, 0, 0));
            y_pred = y_init + delta * (tang_y / dist(tang_x, tang_y, 0, 0));
        }

        if (prints) {
            printf("Prediction: (%+.8le, %+.8le)\n", x_pred, y_pred);
        }

        /* Correction */
        double x_cor, y_cor;
        double a, b, c, d;
        double iteration_distance;
        int counter_cor = 0;
        do {
            a = dx_h(x_pred, y_pred, niub);
            b = dy_h(x_pred, y_pred, niub);
            c = 2 * (x_pred - x_init);
            d = 2 * (y_pred - y_init);

            x_cor = x_pred - (1.0 / (a * d - b * c)) * (d * h(x_pred, y_pred, x0, y0, niub) -
                                                b * ((x_pred - x_init) * (x_pred - x_init) + (y_pred - y_init) * (y_pred - y_init) - delta * delta));
            y_cor = y_pred - (1.0 / (a * d - b * c)) * (a * ((x_pred - x_init) * (x_pred - x_init) + (y_pred - y_init) * (y_pred - y_init) - delta * delta) -
                                                c * h(x_pred, y_pred, x0, y0, niub));

            iteration_distance = dist(x_pred, y_pred, x_cor, y_cor);
            if (prints) {
                printf("Correction: (%+.8le, %+.8le) with distance %+.8le\n", x_cor, y_cor,
                       dist(x_pred, y_pred, x_cor, y_cor));
            }
            x_pred = x_cor;
            y_pred = y_cor;
            counter_cor++;

        } while (iteration_distance >= pre && counter_cor < iter);

        fprintf(fptr, "%6d\t%+.8le\t%+.8le\n", iter_counter, x_pred, y_pred);

        x_init = x_pred;
        y_init = y_pred;
        iter_counter++;

    } while ((dist(x_pred, y_pred, x0, y0) >= delta  || iter_counter < 2) && iter_counter < num);

}

/* ----------------- Auxiliary Functions ------------------ */

double f(double x, double y, const int niub[]) {
    int i;
    double prod = 1;
    for (i = 0; i <= 7; i += 2) {
        prod *= dist(x, y, niub[i], niub[i+1]);
    }
    return prod;
}

double dx_f(double x, double y, const int niub[]) {
    double a = dist(x, y, niub[0], niub[1]);
    double dx_a = dx_dist(x, y, niub[0], niub[1]);
    double b = dist(x, y, niub[2], niub[3]);
    double dx_b = dx_dist(x, y, niub[2], niub[3]);
    double c = dist(x, y, niub[4], niub[5]);
    double dx_c = dx_dist(x, y, niub[4], niub[5]);
    double d = dist(x, y, niub[6], niub[7]);
    double dx_d = dx_dist(x, y, niub[6], niub[7]);
    return dx_a * b * c * d + a * (dx_b * (c * d) + b * (dx_c * d + c * dx_d));
}

double dy_f(double x, double y, const int niub[]) {
    double a = dist(x, y, niub[0], niub[1]);
    double dy_a = dy_dist(x, y, niub[0], niub[1]);
    double b = dist(x, y, niub[2], niub[3]);
    double dy_b = dy_dist(x, y, niub[2], niub[3]);
    double c = dist(x, y, niub[4], niub[5]);
    double dy_c = dy_dist(x, y, niub[4], niub[5]);
    double d = dist(x, y, niub[6], niub[7]);
    double dy_d = dy_dist(x, y, niub[6], niub[7]);
    return dy_a * b * c * d + a * (dy_b * (c * d) + b * (dy_c * d + c * dy_d));
}

double h(double x, double y, double x0, double y0, const int niub[]) {
    return f(x, y, niub) - f(x0, y0, niub);
}

double dx_h(double x, double y, const int niub[]) {
    return dx_f(x, y, niub);
}

double dy_h(double x, double y, const int niub[]) {
    return dy_f(x, y, niub);
}

double dist(double x, double y, double p, double q) {
    return sqrt((x - p) * (x - p) + (y - q) * (y - q));
}

double dx_dist(double x, double y, double p, double q) {
    return (x - p) / dist(x, y, p, q);
}

double dy_dist(double x, double y, double p, double q) {
    return (y - q) / dist(x, y, p, q);
}

/* ----------------- Component Computation ------------------ */

double compute_b(const int niub[], int comp) {

    int i;
    double b = 0;

    if (comp < 0 || comp > 1) {
        printf("Comp (component) must be either 0 or 1!\n\n Exiting...\n");
        exit(2);
    }

    for (i = comp; i <= 7; i += 2) {
        b += niub[i];
    }

    return 0.25 * b;
}

double compute_initial_points(const int niub[], const double distances[], double const b[],
                              int comp, int point, int dist_idx) {

    if (comp < 0 || comp > 1) {
        printf("'comp' (component) must be either 0 or 1!\n\n Exiting...\n");
        exit(2);
    }

    if (point < 0 || point > 3) {
        printf("'point' must be an integer between 0 or 3!\n\n Exiting...\n");
        exit(2);
    }

    if (dist_idx < 0 || dist_idx > 3) {
        printf("'dist_idx' (distances index) must be an integer between 0 or 3!\n\n Exiting...\n");
        exit(2);
    }

    return niub[2*point + comp] + distances[dist_idx] * (niub[2*point + comp] - b[comp]);
}

#include <stdio.h>
#include <math.h>

#define c 0.05
#define rho 2.5

// Define the pendulum ODE as a first-order system
void pendulum_ode(double t, double y[2], double result[2]) {
    double theta = y[0];
    double omega = y[1];
    
    result[0] = omega;
    result[1] = -c * omega - sin(theta) + rho * cos(t);
}

double map_pi(double theta) {
    // Map theta to [-pi, pi)
    double x = fmod(theta, 2.0*M_PI); 
    if (x >= M_PI) {
        x -= 2.0*M_PI;
    }
    return x;
}

void pendulum_orbit(double theta0, double theta_dot0, double* restrict theta, double* restrict theta_dot, int num_points, double dt) {
    double t_max = ((double)num_points) * 2.0*M_PI;
    long n = (long)(t_max / dt) + 1;

    // Solve the ODE using Runge-Kutta 4th order method
    double y[2] = {theta0, theta_dot0};
    double y_temp[2];
    double k1[2], k2[2], k3[2], k4[2];
    double t = 0.0;
    double tcounter = 0.0;
    int k = 0;

    for (long i = 0; i < n; i++) {
        pendulum_ode(t, y, k1);

        for (int j = 0; j < 2; j++) {
            y_temp[j] = y[j] + 0.5 * dt * k1[j];
        }
        pendulum_ode(t + 0.5 * dt, y_temp, k2);

        for (int j = 0; j < 2; j++) {
            y_temp[j] = y[j] + 0.5 * dt * k2[j];
        }
        pendulum_ode(t + 0.5 * dt, y_temp, k3);

        for (int j = 0; j < 2; j++) {
            y_temp[j] = y[j] + dt * k3[j];
        }
        pendulum_ode(t + dt, y_temp, k4);

        for (int j = 0; j < 2; j++) {
            y[j] = y[j] + (dt / 6.0) * (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]);
        }

        t += dt;
        tcounter += dt;
        if (tcounter >= 2.0*M_PI) {
            tcounter -= 2.0*M_PI;
            theta[k] = y[0]; // Map theta to [-pi, pi)
            theta_dot[k] = y[1];
            k += 1;
            if (k == num_points) {
                break;
            }
        }
    }
}
// Runge-Kutta 4th order method 
#include <stdio.h> 
#include <stdlib.h> 
 
// Function to calculate dy/dx = f(x, y) 
double func(double x, double y) { 
    return (x + y); // Example differential equation: dy/dx = x + y 
} 
 
int main() { 
    double x_0, y_0, h, x_n; 
 
    // Input initial conditions and step size 
    printf("Enter the initial value of x and y: \n"); 
    scanf("%lf %lf", &x_0, &y_0); 
 
    // Input step size 
    printf("Enter the step size: "); 
    scanf("%lf", &h); 
 
    // Input the last value of x for computation 
    printf("Enter the last limit: "); 
    scanf("%lf", &x_n); 
 
    // Calculate the number of steps required based on the range [x_0, x_n] 
    int n = (int)((x_n - x_0) / h) + 1; 
 
    // Dynamically allocate memory for the solution table (2 rows: x and y) 
    double **table = (double **)malloc(2 * sizeof(double *)); 
    for (int i = 0; i < 2; i++) { 
        table[i] = (double *)malloc(n * sizeof(double)); 
    } 
 
    // Initialize the first row with the initial conditions 
    table[0][0] = x_0; // First value of x 
    table[1][0] = y_0; // First value of y 
 
    // Apply the Runge-Kutta method to calculate y for each step 
    for (int i = 1; i < n; i++) { 
        double x = table[0][i - 1]; // Current x value 
        double y = table[1][i - 1]; // Current y value 
 
        // Compute the next x value 
        table[0][i] = x + h; 
 
        // Runge-Kutta calculations: 
        double k1 = h * func(x, y); 
        double k2 = h * func(x + h / 2, y + k1 / 2); 
        double k3 = h * func(x + h / 2, y + k2 / 2); 
        double k4 = h * func(x + h, y + k3); 
 
        // Calculate the increment in y using the Runge-Kutta formula 
        double del_y = (k1 + 2 * k2 + 2 * k3 + k4) / 6; 
 
        // Update the next y value 
        table[1][i] = y + del_y; 
    } 
 
    // Print the results in tabular format 
    printf("The solution is given by:\n"); 
    printf(" x           y\n"); 
    for (int i = 0; i < n; i++) { 
        printf("%lf   %lf\n", table[0][i], table[1][i]); 
    } 
 
    // Free the allocated memory to avoid memory leaks 
    for (int i = 0; i < 2; i++) { 
        free(table[i]); 
    } 
    free(table); 
 
    return 0; 
}
 
// Euler’s method 
#include <stdio.h> 
#include <stdlib.h> 
 
// Function defining the differential equation dy/dx = f(x, y) 
double func(double x, double y) { 
    return (x + y); // Example equation: dy/dx = x + y 
} 
 
int main() { 
    double x_0, y_0, h, x_n; 
 
    // Prompt the user for initial conditions (x_0, y_0) 
    printf("Enter the initial value of x and y: \n"); 
    scanf("%lf %lf", &x_0, &y_0); 
 
    // Prompt the user for the step size (h) 
    printf("Enter the step size: "); 
    scanf("%lf", &h); 
 
    // Prompt the user for the last x value (x_n) up to which to solve the equation 
    printf("Enter the last limit: "); 
    scanf("%lf", &x_n); 
 
    // Calculate the number of steps needed for the given range [x_0, x_n] 
    int n = (int)((x_n - x_0) / h) + 1; 
 
    // Dynamically allocate memory for a 2D table (2 rows: x and y values) 
    double **table = (double **)malloc(2 * sizeof(double *)); 
    for (int i = 0; i < 2; i++) { 
        table[i] = (double *)malloc(n * sizeof(double)); 
    } 
 
    // Initialize the table with the initial conditions 
    table[0][0] = x_0; // Set the initial x value 
    table[1][0] = y_0; // Set the initial y value 
 
    // Perform the numerical solution using Euler's method 
    for (int i = 1; i < n; i++) { 
        double x = table[0][i - 1]; // Current x value 
        double y = table[1][i - 1]; // Current y value 
        table[0][i] = x + h; // Compute the next x value 
        table[1][i] = y + h * func(x, y); // Compute the next y value using Euler's method 
    } 
 
    // Print the results in tabular format 
    printf("The solution is given by:\n"); 
    printf(" x           y\n"); 
    for (int i = 0; i < n; i++) { 
        printf("%lf   %lf\n", table[0][i], table[1][i]); 
    } 
 
    // Free the dynamically allocated memory to avoid memory leaks 
    for (int i = 0; i < 2; i++) { 
        free(table[i]); // Free each row 
    } 
    free(table); // Free the top-level pointer 
 
    return 0; 
}
// Simpson’s 3/8th rule 
#include <stdio.h> 
#include <conio.h> 
#include <math.h> 
 
// Function to evaluate the integrand 
double function(double x) { 
    return 1 / ((x * x) + 1); // The function to integrate 
} 
 
int main() { 
    int n; 
    // Input the number of terms 
    printf("Enter the number of terms: "); 
    scanf("%d", &n); // Number of terms or subintervals for Simpson's 3/8th rule 
 
    int upperLimit, lowerLimit; 
 
    // Input the integration limits 
    printf("Enter the lower limit and upper limit: "); 
    scanf("%d%d", &lowerLimit, &upperLimit); 
 
    // Calculate the width of each subinterval 
    double h = (double)(upperLimit - lowerLimit) / (n - 1); 
 
    double y[n]; // Array to store function values at each subinterval point 
 
    // Evaluate the function at each subinterval 
    for (int i = 0; i < n; i++) { 
        double x = (h * i) + lowerLimit; // x-coordinate for each point 
        y[i] = function(x); 
        // Store function value in y[i] 
        printf("x_%d : %lf    y_%d : %lf\n", i, x, i, y[i]); 
    } 
 
    // Initialize the result with the first and last function values 
    double ans = y[0] + y[n - 1]; 
 
    // Apply Simpson's 3/8th rule 
    for (int i = 1; i < n - 1; i++) { 
        if (i % 3 == 0) { 
            ans += 2 * y[i]; 
        } else { 
            ans += 3 * y[i]; 
        } 
    } 
 
    // Finalize the result 
    ans *= (3 * h / 8); 
 
    // Print the result 
    printf("\nThe value of the integral is %lf", ans); 
 
    return 0; 
}
// Simpson’s 1/3rd rule 
#include <stdio.h> 
#include <conio.h> 
#include <math.h> 
 
// Function to evaluate the integrand 
double function(double x) { 
    return 1 / ((x * x) + 1); // The function to integrate 
} 
 
int main() { 
    int n; 
    // Input the number of terms 
    printf("Enter the number of terms: "); 
    scanf("%d", &n); // Number of terms or subintervals for Simpson's 1/3rd rule 
 
    int upperLimit, lowerLimit; 
 
    // Input the integration limits 
    printf("Enter the lower limit and upper limit: "); 
    scanf("%d%d", &lowerLimit, &upperLimit); 
 
    // Calculate the width of each subinterval 
    double h = (double)(upperLimit - lowerLimit) / (n - 1); 
 
    double y[n]; // Array to store function values at each subinterval point 
 
    // Evaluate the function at each subinterval 
    for (int i = 0; i < n; i++) { 
        double x = (h * i) + lowerLimit; // x-coordinate for each point 
        y[i] = function(x); 
        // Store function value in y[i] 
        printf("x_%d : %lf    y_%d : %lf\n", i, x, i, y[i]); 
    } 
 
    // Initialize the result with the first and last function values 
    double ans = y[0] + y[n - 1]; 
 
    // Apply Simpson's 1/3rd rule 
    for (int i = 1; i < n - 1; i++) { 
        if (i % 2 == 0) { 
            ans += 2 * y[i]; 
        } else { 
            ans += 4 * y[i]; 
        } 
    } 
 
    // Finalize the result 
    ans *= (h / 3); 
 
    // Print the result 
    printf("\nThe value of the integral is %lf", ans); 
 
    return 0; 
}
// Trapezoidal rule 
#include <stdio.h> 
#include <conio.h> 
#include <math.h> 
 
// Function to evaluate the integrand 
double function(double x) { 
    return 1 / ((x * x * x) + 10); // Function to integrate 
} 
 
int main() { 
    int n; 
    printf("Enter the number of terms: "); 
    scanf("%d", &n); // Number of terms/subintervals for the trapezoidal rule 
 
    int upperLimit, lowerLimit; 
    printf("Enter the lower and upper limits: "); 
    scanf("%d%d", &lowerLimit, &upperLimit); // Integration limits 
 
    // Calculate the width of each subinterval 
    double h = (double)(upperLimit - lowerLimit) / (n - 1); 
 
    double y[n]; // Array to store function values 
 
    // Evaluate the function at each subinterval point 
    for (int i = 0; i < n; i++) { 
        double x = (h * i) + lowerLimit; // x-coordinate for each point 
        y[i] = function(x); 
        // Store function value in y[i] 
        printf("x_%d : %lf    y_%d : %lf\n", i, x, i, y[i]); 
    } 
 
    // Initialize the result with the first and last function values 
    double ans = y[0] + y[n - 1]; 
 
    // Apply the trapezoidal rule for inner points 
    for (int i = 1; i < n - 1; i++) { 
        ans += (2 * y[i]); 
    } 
 
    // Finalize the result 
    ans *= (h / 2); 
 
    // Print the result 
    printf("\n\nThe value of the integral is %lf", ans); 
 
    return 0; 
}
 
// Langrange’s interpolation method 
#include <stdio.h> 
#include <math.h> 
 
int main() { 
    int n; 
    // Input the number of data points 
    printf("Enter the number of entries: "); 
    scanf("%d", &n); 
 
    double x_val;       // x-value to find the corresponding y 
    double y_val = 0;   // Initialize y_val to store the interpolated result 
 
    int x[n], y[n]; 
 
    // Input x and y values 
    printf("Enter the values of x and y:\n"); 
    for (int i = 0; i < n; i++) { 
        scanf("%d", &x[i]); 
        scanf("%d", &y[i]); 
    } 
 
    // Input the x-value for interpolation 
    printf("Enter the value of x to find y: "); 
    scanf("%lf", &x_val); 
 
    // Perform Lagrange Interpolation 
    for (int i = 0; i < n; i++) { 
        double multi = 1.0; // Numerator product term 
        double div = 1.0;   // Denominator product term 
 
        for (int j = 0; j < n; j++) { 
            if (j != i) { // Exclude the current point from calculation 
                multi *= (x_val - x[j]); // Numerator 
                div *= (x[i] - x[j]);   // Denominator 
            } 
        } 
 
        // Add the contribution of the i-th term to y_val 
        y_val += (multi / div) * y[i]; 
    } 
 
    // Output the interpolated y-value 
    printf("y(%lf) = %lf", x_val, y_val); 
 
    return 0; 
}
// Newton Backward Interpolation method 
#include <stdio.h> 
#include <math.h> 
 
int factorial(int a) { 
    int fact = 1; 
    while (a > 1) { 
        fact *= a; 
        a--; 
    } 
    return fact; 
} 
 
int main() { 
    int n; 
 
    printf("Enter the number of values: "); 
    scanf("%d", &n); 
 
    int x[n], y[n]; 
    int polation[n]; 
    int y_table[n][n]; 
 
    // Input the values of x and y table 
    printf("Enter the values of x and y table:\n"); 
    for (int i = 0; i < n; i++) { 
        scanf("%d%d", &x[i], &y[i]); 
    } 
 
    int x_val; 
 
    // Input the x value for interpolation 
    printf("Enter the value whose y to be found: "); 
    scanf("%d", &x_val); 
 
    double h = x[1] - x[0];  // Calculate step size (h) 
    double u = (double)(x_val - x[n-1]) / h;  // Calculate (x - xn)/h for backward interpolation 
 
    polation[0] = y[n-1];  // Initialize interpolation with the last value of y 
 
    // Fill in the y_table with initial y values 
    for (int i = 0; i < n; i++) { 
        y_table[i][0] = y[i]; 
    } 
 
    // Calculate the backward difference table 
    for (int i = 1; i < n; i++) { 
        for (int j = i; j < n; j++) { 
            y_table[j][i] = y_table[j][i - 1] - y_table[j - 1][i - 1]; 
        } 
        polation[i] = y_table[n - 1][i];  // Store backward differences for interpolation 
    } 
 
    // Display the backward difference table 
    printf("The backward difference table is:\n"); 
    for (int i = 0; i < n; i++) { 
        for (int j = 0; j <= i; j++) { 
            printf("%d ", y_table[i][j]); 
        } 
        printf("\n"); 
    } 
 
    // Perform interpolation using Newton's backward interpolation formula 
    double ans = polation[0];  // Initialize with y[n-1] 
    double multiply = 1.0;  // Variable to accumulate product terms 
 
    for (int i = 1; i < n; i++) { 
        multiply *= (u + (i - 1));  // Calculate product (u)(u+1)...(u+i-1) 
        ans += (multiply * polation[i]) / factorial(i);  // Add ith term of interpolation formula 
    } 
 
    // Print the interpolated result 
    printf("f(%d) = %.3lf\n", x_val, ans); 
 
    return 0; 
}
// Newton Forward Interpolation method 
#include <stdio.h> 
#include <math.h> 
 
int factorial(int a) { 
    int fact = 1; 
    while (a > 1) { 
        fact *= a; 
        a--; 
    } 
    return fact; 
} 
 
int main() { 
    int n; 
    printf("Enter the number of values: "); 
    scanf("%d", &n); 
 
    int x[n], y[n]; 
    int polation[n]; 
    int y_table[n][n]; 
 
    printf("Enter the values of x and y table:\n"); 
    for (int i = 0; i < n; i++) { 
        scanf("%d%d", &x[i], &y[i]); 
    } 
 
    int x_val; 
     
    printf("Enter the value whose y is to be found: "); 
    scanf("%d", &x_val); 
 
    double h = x[1] - x[0];  // Calculate step size (h) 
    double u = (double)(x_val - x[0]) / h;  // Calculate (x - x0)/h 
 
    polation[0] = y[0];  // Initialize interpolation polynomial 
 
    // Fill in the y_table with initial y values 
    for (int i = 0; i < n; i++) { 
        y_table[i][0] = y[i]; 
    } 
 
    // Calculate the forward difference table 
    for (int i = 1; i < n; i++) { 
        for (int j = 0; j < n - i; j++) { 
            y_table[j][i] = y_table[j + 1][i - 1] - y_table[j][i - 1]; 
        } 
        polation[i] = y_table[0][i];  // Store forward differences for interpolation 
    } 
 
    // Display the forward difference table 
    printf("The forward difference table is:\n"); 
    for (int i = 0; i < n; i++) { 
        for (int j = 0; j < n - i; j++) { 
            printf("%d ", y_table[i][j]); 
        } 
        printf("\n"); 
    } 
 
    // Perform interpolation using Newton's forward interpolation formula 
    double ans = polation[0];  // Initial y0 value 
    double multiply = 1.0;  // Variable to accumulate product terms 
 
    for (int i = 1; i < n; i++) { 
        multiply *= (u - (i - 1));  // Calculate product (u)(u-1)...(u-i+1) 
        ans += (multiply * polation[i]) / factorial(i);  // Add ith term of interpolation formula 
    } 
 
    // Print the interpolated result 
    printf("f(%d) = %.3lf\n", x_val, ans); 
 
    return 0; 
}
// Gauss Jordan
#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // Include math.h for fabs()
 
int main() {
    int n;
    printf("Enter the number of equations: ");
    scanf("%d", &n);
 
    float a[n][n + 1]; // Augmented matrix
    float ans[n];       // Array to store the solutions
    float ratio;        // Ratio used for eliminating variables
 
    // Input the coefficients of the variables
    printf("Enter the coefficients of the variables:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            scanf("%f", &a[i][j]);
        }
    }
 
    // Display the augmented matrix
    printf("The augmented matrix is:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            printf("%0.3f ", a[i][j]);
        }
        printf("\n");
    }
 
    // Gaussian Elimination with Row Swapping for stability
    for (int i = 0; i < n; i++) {
        // Pivoting: Find the row with the maximum element in the current column
        if (a[i][i] == 0.0) {
            int max_row = i;
            // Find the row with the maximum element in the current column
            for (int k = i + 1; k < n; k++) {
                if (fabs(a[k][i]) > fabs(a[max_row][i])) {
                    max_row = k;
                }
            }
 
            // Swap the current row with the row having the maximum element in the column
            if (max_row != i) {
                for (int j = 0; j <= n; j++) {
                    float temp = a[i][j];
                    a[i][j] = a[max_row][j];
                    a[max_row][j] = temp;
                }
            }
        }
 
        // If the pivot element is zero after row swapping, exit the program
        if (a[i][i] == 0.0) {
            printf("The system has no unique solution.\n");
            exit(0);
        }
 
        // Perform elimination to create an upper triangular matrix
        for (int j = i + 1; j < n; j++) {
            ratio = a[j][i] / a[i][i];
            for (int k = 0; k <= n; k++) {
                a[j][k] = a[j][k] - ratio * a[i][k];
            }
        }
    }
 
    // Display the upper triangular matrix
    printf("The augmented upper triangular matrix is:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            printf("%0.3f ", a[i][j]);
        }
        printf("\n");
    }
 
    // Back substitution to find the solutions
    ans[n - 1] = a[n - 1][n] / a[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--) {
        ans[i] = a[i][n];
        for (int j = n - 1; j > i; j--) {
            ans[i] -= a[i][j] * ans[j];
        }
        ans[i] = ans[i] / a[i][i];
    }
 
    // Display the solution
    printf("The solution is:\n");
    for (int i = 0; i < n; i++) {
        printf("x%d = %0.3f\n", i + 1, ans[i]);
    }
 
    return 0;
}
// Gauss elinination
#include <stdio.h>
#include <stdlib.h>
 
int main() {
    int n;
    printf("Enter the number of equations: ");
    scanf("%d", &n);
 
    float a[n][n + 1]; // Augmented matrix
    float ans[n];       // Array to store the solutions
    float ratio;        // Ratio used for eliminating variables
 
    // Input the coefficients of the variables
    printf("Enter the coefficients of the variables:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            scanf("%f", &a[i][j]);
        }
    }
 
    // Display the augmented matrix
    printf("The augmented matrix is:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            printf("%0.3f ", a[i][j]);
        }
        printf("\n");
    }
 
    // Gaussian Elimination with Row Swapping
    for (int i = 0; i < n; i++) {
        // Pivoting: Find the row with the maximum element in the current column
        if (a[i][i] == 0.0) {
            int max_row = i;
            for (int k = i + 1; k < n; k++) {
                if (fabs(a[k][i]) > fabs(a[max_row][i])) {
                    max_row = k;
                }
            }
 
            // Swap the current row with the row having the maximum element
            if (max_row != i) {
                for (int j = 0; j <= n; j++) {
                    float temp = a[i][j];
                    a[i][j] = a[max_row][j];
                    a[max_row][j] = temp;
                }
            }
        }
 
        // If the pivot element is zero after row swapping, exit the program
        if (a[i][i] == 0.0) {
            printf("The system has no unique solution.\n");
            exit(0);
        }
 
        // Perform elimination to create upper triangular matrix
        for (int j = i + 1; j < n; j++) {
            ratio = a[j][i] / a[i][i];
            for (int k = 0; k <= n; k++) {
                a[j][k] = a[j][k] - ratio * a[i][k];
            }
        }
    }
 
    // Display the upper triangular matrix
    printf("The augmented upper triangular matrix is:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            printf("%0.3f ", a[i][j]);
        }
        printf("\n");
    }
 
    // Back substitution to find the solutions
    ans[n - 1] = a[n - 1][n] / a[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--) {
        ans[i] = a[i][n];
        for (int j = n - 1; j > i; j--) {
            ans[i] -= a[i][j] * ans[j];
        }
        ans[i] = ans[i] / a[i][i];
    }
 
    // Display the solution
    printf("The solution is:\n");
    for (int i = 0; i < n; i++) {
        printf("x%d = %0.3f\n", i + 1, ans[i]);
    }
 
    return 0;
}
// Newton Raphson Method 
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
 
#define f(x) (3 * (x) - cos(x) - 1) 
#define g(x) (3 + sin(x)) 
 
int main() { 
    float x0_1, x0_2, x0, x1, f0, f1, g0, e;     
    int i = 1, N; 
 
    printf("Enter the first guess: "); 
    scanf("%f", &x0_1); 
 
    printf("Enter the second guess: "); 
    scanf("%f", &x0_2); 
 
    printf("Enter the error value: "); 
    scanf("%f", &e); 
 
    printf("Enter the maximum iterations: "); 
    scanf("%d", &N); 
 
    if (fabs(f(x0_1)) < fabs(f(x0_2))) { 
        x0 = x0_1; 
    } else { 
        x0 = x0_2; 
    } 
 
    printf("Initial guess is: %f\n", x0); 
    printf("\nStep\t\tx0\t\tf(x0)\t\tx1\t\tf(x1)\n"); 
 
    do { 
        g0 = g(x0); 
        f0 = f(x0); 
 
        if (g0 == 0.0) { 
            printf("Mathematical Error: Derivative is zero.\n"); 
            exit(0); 
        } 
 
        x1 = x0 - f0 / g0; 
        f1 = f(x1); 
 
        printf("%d\t\t%.6f\t\t%.6f\t\t%.6f\t\t%.6f\n", i, x0, f0, x1, f1); 
 
        x0 = x1; 
        i++; 
 
        if (i > N) { 
            printf("Not Convergent: Exceeded maximum iterations.\n"); 
            exit(0); 
        } 
 
    } while (fabs(f1) > e); 
 
    printf("\nApproximate root is: %f\n", x1); 
 
    return 0; 
}
// Secant Method 
#include <stdio.h> 
#include <math.h> 
#include <stdlib.h> 
 
#define f(x) (x * x * x - 2 * x - 5) 
 
int main() { 
    float x0, x1, x2, f0, f1, f2, e;     
    int i = 1, N; 
 
    printf("Enter the first guess: "); 
    scanf("%f", &x0); 
 
    printf("Enter the second guess: "); 
    scanf("%f", &x1); 
     
    f0 = f(x0); 
    f1 = f(x1); 
 
    printf("Enter the error value: "); 
    scanf("%f", &e); 
 
    printf("Enter the maximum iterations: "); 
    scanf("%d", &N); 
 
    printf("\nStep\tx0\t\tx1\t\tx2\n"); 
 
    do { 
        f0 = f(x0); 
        f1 = f(x1); 
 
        if (fabs(f1 - f0) < 1e-10) { 
            printf("Mathematical error: f(x0) and f(x1) are too close.\n"); 
            exit(0); 
        } 
 
        x2 = x1 - (x1 - x0) * f1 / (f1 - f0); 
        f2 = f(x2); 
 
        printf("%d\t%.6f\t\t%.6f\t\t%.6f\n", i, x0, x1, x2); 
 
        x0 = x1; 
        f0 = f1; 
        x1 = x2; 
        f1 = f2; 
 
        i++; 
 
        if (fabs(x1 - x0) < e) { 
            break; 
        } 
 
        if (i > N) { 
            printf("Not convergent after %d iterations.\n", N); 
            exit(0); 
        } 
 
    } while (1); 
 
    printf("Approximate root is: %.6f\n", x2); 
 
    return 0; 
}
// Regula Falsi Method 
#include <stdio.h> 
#include <math.h> 
 
float f(float x) { 
     return cos(x) - x * exp(x); 
} 
int main() { 
    float x0, x1, x2, f0, f1, f2, e;     
    int i = 1; 
 
    printf("Enter the first guess: "); 
    scanf("%f", &x0); 
    printf("Enter the second guess: "); 
    scanf("%f", &x1); 
    printf("Enter the error value: "); 
    scanf("%f", &e); 
 
    f0 = f(x0);     
    f1 = f(x1); 
 
    if (f0* f1 > 0) { 
        printf("Given guesses are incorrect. The function does not change sign over the interval.\n"); 
        return 1; 
    } 
 
    printf("\nStep\t\tx0\t\tx1\t\tx2\n"); 
    do { 
        x2 = x0 - (f1 * (x0 - x1)) / (f0 – f1); 
        printf("%d\t\t%.4f\t\t%.4f\t\t%.4f\n", i, x0, x1, x2); 
        f2= f(x2); 
 
        if (f1 * f2 < 0) { 
            x1 = x2; 
            f1 = f2; 
        } else { 
            x0 = x2; 
            f0 = f2; 
        } 
        i++; 
    } while (fabs(f2) > e); 
 
    printf("Total number of iterations: %d\n", i - 1); 
    printf("Approximate root: %.4f\n", x2); 
    return 0; 
}
// Bisection Method 
#include <stdio.h> 
#include <math.h> 
 
float f(float x) { 
    return x * log10(x) - 1.2; 
} 
 
int main() {     
    float x0, x1, x2, e;     
    int i = 1; 
 
    printf("Enter the first guess: ");     
    scanf("%f", &x0); 
 
    printf("Enter the second guess: ");     
    scanf("%f", &x1); 
 
    printf("Enter the error value: ");     
    scanf("%f", &e); 
 
    if (f(x0) * f(x1) > 0) { 
        printf("The function does not change sign over the interval [%.2f, %.2f].\n", x0, x1); 
        return 1; 
    } 
 
    printf("\nStep\t\tx0\t\tx1\t\tx2\t\terror\n"); 
     
    do { 
        x2 = (x0 + x1) / 2; 
        printf("%d\t\t%.4f\t%.4f\t%.4f\t%.4f\n", i, x0, x1, x2, fabs(x1 - x0)); 
 
        if (f(x0) * f(x2) < 0) { 
            x1 = x2; 
        } else if (f(x0) * f(x2) > 0) { 
            x0 = x2; 
        } 
        i++;   
    } while (fabs(x1 - x0) > e); 
 
    printf("Approximate root is: %.4f\n", x2); 
 
    return 0; 
}

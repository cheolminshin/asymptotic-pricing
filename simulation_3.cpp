#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>
#include <chrono>

using namespace std ;
using namespace std::chrono ;

double nsamples = 1e3 ;
double ntimes = 1e3 ;

//Normal distribution sample generator using Box-Muller transform
double normal()
{
    static int count = 0 ;
    static double storage ;
    if(count==0){
        double u1, u2, z1, z2 ;
        do{
            u1 = rand()/(double)RAND_MAX ;
            u2 = rand()/(double)RAND_MAX ;
            z1 = sqrt(-2*log(u1))*cos(2*M_PI*u2) ;
            z2 = sqrt(-2*log(u1))*sin(2*M_PI*u2) ;
        }while(u1==0);
        storage = z2 ;
        count = 1 ;
        return z1 ;
    }
    else{
        count = 0 ;
        return storage ;
    }
}

// Heston function

double asian_heston_original(double nsamples,
                            double ntimes,
                            double tau,
                            double s, 
                            double nu, 
                            double kappa, 
                            double theta, 
                            double xi,
                            double rho,
                            double r,
                            double K)
{
    double dw1, dw2 ;               // normal distribution random variable container
    double integs = 0.0 ;
    double out = 0.0 ;              // output
    double h = tau / ntimes ;       // time increment
    double sqrth = sqrt(h) ;        // square root of h
    double rho2 = sqrt(1-pow(rho,2)) ;
    double st, nut ;                // current value container of stock price and volatility

    int nonzeropayoff = 0 ;

    for (int i=0; i<nsamples; i++)
    {
        st = s ;                    // initialize the current stock price
        nut= nu ;                   // initialize the current volatility
        integs = 0.0 ;              // initialize the time integral of S
        for (int j=0; j<ntimes; j++)
        {
            dw1 = sqrth * normal() ;    // assign Brownian time increment of h
            dw2 = sqrth * normal() ;    // assign Brownian time increment of h

            st *= exp((r - 0.5 * nut) * h + sqrt(nut) * dw1);
            nut += kappa*(theta-nut)*h + xi*sqrt(nut)*(rho*dw1 + rho2*dw2) ;      // Update the current volatility by Euler-Maruyama
            nut = max(nut, 1e-8) ;

            integs += st ;      // Update the integral value
        }

        integs /= ntimes ;           // Adjust the time scale
        integs = max(integs-K, 0.0) ;
        out += integs ;         // Update the sample value
    }
    out *= exp(-r*tau)/nsamples ;           // result.
    return out ;
}

// page 10 of 16, proposition 4, sigma_A bar^2 (tau)
double asian_volatility_heston(double tau, 
                                double x,
                                double nu, 
                                double kappa,
                                double theta,
                                double xi,
                                double rho)
{
    double A = nu/3.0 ;
    double B = rho*xi/6.0 ;
    double C = pow(rho*xi,2)/(48.0 * nu) + (2.0-5.0*pow(rho,2))*pow(xi,2) /(72.0*nu) ;
    double D = (rho*xi*nu)/24.0 + kappa*(theta-nu)/12.0 + (pow(rho,2)-4)*pow(xi,2)/288.0 ;
    double E = pow(rho*xi,2)/144.0 
                    + kappa*rho*xi*(3.0*theta-2.0*nu)/(144.0*nu)
                    + rho*pow(xi,2)*(pow(rho,2)*xi-4.0*xi-9*pow(rho,2)+12.0)/(1152.0*nu) ;
    
    double out = A + B*x + C*pow(x,2) + D*tau + E*x*tau ;
    return out ;
}

double asian_heston_second(double nsamples,
                            double tau,
                            double s, 
                            double nu, 
                            double kappa, 
                            double theta, 
                            double xi,
                            double rho,
                            double r,
                            double K)
{
    double out = 0.0 ;
    double x = log(K/s) - r*tau ;       // x = log(K/s exp(r tau))
    double sqrttau = sqrt(tau) ;
    double sigma = sqrt(asian_volatility_heston(tau,x,nu,kappa,theta,xi,rho)) ;
    double w ;
    for (int i=0; i<nsamples; i++)
    {
        w = normal() ;
        out += max(1.0+ sigma* sqrttau* w -K/s, 0.0) ;
    }
    out *= (s/nsamples) ;
    return out ;
}


// SABR

double asian_sabr_original(double nsamples,
                            double ntimes,
                            double tau,
                            double s, 
                            double nu, 
                            double alpha, 
                            double beta,
                            double rho,
                            double r,
                            double K)
{
    double dw1, dw2 ;               // normal distribution random variable container
    double integs = 0.0 ;
    double out = 0.0 ;              // output
    double h = tau / ntimes ;       // time increment
    double sqrth = sqrt(h) ;        // square root of h
    double rho2 = sqrt(1-pow(rho,2)) ;
    double ft, st, nut ;                // current value container of forward price, stock price, and volatility
    double t ;                      // current time

    //int nonzeropayoff = 0 ;

    for (int i=0; i<nsamples; i++)
    {
        ft = s * exp(r*tau) ;       // initialize the current forward price
        nut= nu ;                   // initialize the current volatility
        integs = 0.0 ;              // initialize the time integral of S
        t = 0.0 ;
        for (int j=0; j<ntimes; j++)
        {
            t += h ;
            dw1 = sqrth * normal() ;    // assign Brownian time increment of h
            dw2 = sqrth * normal() ;    // assign Brownian time increment of h

            ft += nut * pow(ft,beta) * dw1;                 // Update the current forward price by Euler-Maruyama
            nut += alpha * nut *(rho*dw1 + rho2*dw2) ;      // Update the current volatility by Euler-Maruyama
            ft = max(ft, 1e-8) ;        // prevent being negative
            nut = max(nut, 1e-8) ;      // prevent being negative

            st = ft * exp(-r*(tau-t)) ;                     // Update the current St

            integs += st ;              // Update the integral value
        }

        integs /= ntimes ;              // Adjust the time scale
        integs = max(integs-K, 0.0) ;
        out += integs ;                 // Update the sample value
    }
    out *= exp(-r*tau)/nsamples ;           // result.
    return out ;
}

// page 10 of 16, proposition 4, sigma_A bar^2 (tau)
double asian_volatility_sabr(double tau, 
                                double x,
                                double z,
                                double nu, 
                                double alpha,
                                double beta,
                                double rho)
{
    double z1 = pow(z,beta-1.0) ;
    double z2 = pow(z,2.0*(beta-1.0)) ;
    double z3 = pow(z,3.0*(beta-1.0)) ;
    double z4 = pow(z,4.0*(beta-1.0)) ;

    double A = pow(nu,2.0) * z2 / 3.0 ;
    double B = alpha*rho*nu*z1/3.0 + (beta-1.0)*pow(nu,2.0)*z2/3.0 ;
    double C = (4.0-6.0*pow(beta,2.0)+3.0*pow(rho,2.0))*pow(alpha,2.0)/36.0
                    + (beta-1.0)*alpha*rho*nu*z1/6.0
                    + 5.0*pow(nu*(1.0-beta),2.0)*z2/36.0 ;
    double D = (2.0-3.0*pow(rho,2.0))*pow(alpha*nu,2.0)*z2/72.0
                    + beta*alpha*rho*pow(nu,3.0)*z3/12.0
                    + pow(1.0-beta,2.0)*pow(nu,4.0)*z4/12.0 ;
    double E = (2.0-3.0*pow(rho,2.0))*pow(alpha,3.0)*rho*nu*z1/144.0
                    + (-beta*pow(rho,2.0)+3.0*pow(rho,2)+2.0*beta-2.0)*pow(alpha*nu,2.0)*z2/48.0
                    + (1.0-31.0*beta)*alpha*rho*pow(nu,3.0)*(1.0-beta)*z3/144.0
                    + 7.0*pow(beta-1.0,3.0)*pow(nu,4.0)*z4/144.0 ;
    
    double out = A + B*x + C*pow(x,2) + D*tau + E*x*tau ;
    return out ;
}

double asian_sabr_second(double nsamples,
                            double tau,
                            double z,
                            double nu,
                            double alpha, 
                            double beta,
                            double rho,
                            double r,
                            double K)
{
    double out = 0.0 ;
    double x = log(K/z) ;       // x = log(K/z)
    double sqrttau = sqrt(tau) ;
    double sigma = sqrt(asian_volatility_sabr(tau,x,z,nu,alpha,beta,rho)) ;
    double w ;
    //cout << "Asian Heston second SIGMA: " << sigma << endl ;
    for (int i=0; i<nsamples; i++)
    {
        w = normal() ;
        out += max(z*exp(-r*tau)*(sigma*sqrttau*w+1.0)-K, 0.0) ;
    }
    out /= nsamples ;
    return out ;
}

string format_time(double seconds)
{
    int mins = static_cast<int>(seconds) / 60;
    int secs = static_cast<int>(seconds) % 60;
    return to_string(mins) + " minutes " + to_string(secs) + " seconds";
}

int main()
{
    double table_of_maturity[8] = {1.0/12.0, 1.0/6.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
    double table_of_strike[7] = {85.0, 90.0, 95.5, 100.0, 105.0, 110.0, 115.0} ;
    int num_of_maturity = sizeof(table_of_maturity) / sizeof(table_of_maturity[0]) ;
    int num_of_K = sizeof(table_of_strike) / sizeof(table_of_strike[0]) ;

    // Temporary variables
    double tau, K, s ;
    // Heston parameters
    double kappa1, theta1, xi1, rho1, r1, nu01, s1 ;
    // SABR parameters
    double z2, nu02, alpha2, beta2, rho2, r2 ;

    // Number of models. Here we have 2 method; Monte-Carlo and Modified method.
    int num_of_model = 2 ;


    kappa1 = 1.3 ;
    theta1 = 0.5 ;
    xi1 = 0.48 ;
    rho1 = -0.6 ;
    r1 = 0.01 ;
    nu01 = 0.3 ;
    s1 = 100.0 ;

    z2 = 100.0 ;
    nu02 = 0.6 ;
    alpha2 = 0.2 ;
    beta2 = 0.7 ;
    rho2 = -0.6 ;
    r2 = 0.01 ;

    int total_iterations ;
    auto start_time = steady_clock::now() ;

    // Table 1
    // Heston model, various strike K's and various maturity tau's
    std::cout << "[Heston model simulation started]" << endl ;
    // Memory Allocation
    double*** out_heston = new double**[num_of_model] ;
    for (int i = 0; i < num_of_model; i++)
    {
        out_heston[i] = new double*[num_of_K] ;
        {
            for (int j = 0; j<num_of_K; j++)
            {
                out_heston[i][j] = new double[num_of_maturity] ;
                for (int k = 0 ; k < num_of_maturity; k++)
                {
                    out_heston[i][j][k] = 0.0 ;
                }
            }
        }
    }
    // Computation Start
    start_time = steady_clock::now();
    total_iterations = num_of_K * num_of_maturity ;
    for (int i=0; i<num_of_K; i++)
    {
        K = table_of_strike[i] ;
        for (int j=0; j<num_of_maturity; j++)
        {
            tau = table_of_maturity[j] ;

            auto now = steady_clock::now();
            duration<double> elapsed = now - start_time;
            double progress = ((double)(i*num_of_K) + (double)(j + 1)) / total_iterations;
            double elapsed_seconds = elapsed.count();
            double estimated_total_time = elapsed_seconds / progress;
            double remaining_time = estimated_total_time - elapsed_seconds;

            std::cout << fixed << setprecision(2);
            std::cout << "\r[Progression] " << i*num_of_K + (j + 1) << " / " << total_iterations;
            std::cout << " (" << progress * 100 << "%)";
            std::cout << " | Elapsed: " << format_time(elapsed_seconds);
            std::cout << " | Remaining: " << format_time(remaining_time) ;
            std::cout.flush();

            out_heston[0][i][j] = asian_heston_original(nsamples,ntimes,tau,s1,nu01,kappa1,theta1,xi1,rho1,r1,K) ;
            out_heston[1][i][j] = asian_heston_second(nsamples,tau,s1,nu01,kappa1,theta1,xi1,rho1,r1,K) ;
        }
    }
    std::cout << "\n[Heson model simulation finished]" << endl;
    // File writing-1
    ofstream file_heston("output_table1(heston).csv");
    if (!file_heston.is_open()) {
        cerr << "Unable to open file!" << endl;
        return 1;
    }

    for (int i = 0; i < num_of_maturity; ++i) {
        for (int j = 0; j < num_of_K; ++j) {
            file_heston << out_heston[0][j][i];
            if (j < num_of_K) file_heston << ",";
        }
        file_heston << "\n";
        for (int j = 0; j < num_of_K; ++j) {
            file_heston << out_heston[1][j][i];
            if (j < num_of_K) file_heston << ",";
        }
        file_heston << "\n";
    }
    file_heston.close();
    std::cout << "Simulation result for Heston has been saved as CSV file!" << endl;
    std::cout << endl;

    // Memory Clear
    for (int i = 0; i < num_of_model; ++i) {
        for (int j = 0; j < num_of_K; ++j) {
            delete[] out_heston[i][j];
        }
        delete[] out_heston[i];
    }
    delete[] out_heston;
    // First Experiment is terminated.



    // Table 2
    // SABR model, various strike K's and various maturity tau's
    std::cout << "[SABR model simulation started]" << endl ;
    // Memory Allocation
    double*** out_sabr = new double**[num_of_model] ;
    for (int i = 0; i < num_of_model; i++)
    {
        out_sabr[i] = new double*[num_of_K] ;
        {
            for (int j = 0; j<num_of_K; j++)
            {
                out_sabr[i][j] = new double[num_of_maturity] ;
                for (int k = 0 ; k < num_of_maturity; k++)
                {
                    out_sabr[i][j][k] = 0.0 ;
                }
            }
        }
    }
    // Computation Start
    start_time = steady_clock::now();
    total_iterations = num_of_K * num_of_maturity ;
    for (int i=0; i<num_of_K; i++)
    {
        K = table_of_strike[i] ;
        for (int j=0; j<num_of_maturity; j++)
        {
            tau = table_of_maturity[j] ;
            s = z2 * exp(-r2*tau) ;

            auto now = steady_clock::now();
            duration<double> elapsed = now - start_time;
            double progress = ((double)(i*num_of_K) + (double)(j + 1)) / total_iterations;
            double elapsed_seconds = elapsed.count();
            double estimated_total_time = elapsed_seconds / progress;
            double remaining_time = estimated_total_time - elapsed_seconds;

            std::cout << fixed << setprecision(2);
            std::cout << "\r[Progression] " << i*num_of_K + (j + 1) << " / " << total_iterations;
            std::cout << " (" << progress * 100 << "%)";
            std::cout << " | Elapsed: " << format_time(elapsed_seconds);
            std::cout << " | Remaining: " << format_time(remaining_time) ;
            std::cout.flush();

            out_sabr[0][i][j] = asian_sabr_original(nsamples,ntimes,tau,s,nu02,alpha2,beta2,rho2,r2,K) ;
            out_sabr[1][i][j] = asian_sabr_second(nsamples,tau,s,nu02,alpha2,beta2,rho2,r2,K) ;
        }
    }
    std::cout << "\n[SABR model simulation finished]" << endl;
    // File writing-1
    ofstream file_sabr("output_table2(SABR).csv");
    if (!file_sabr.is_open()) {
        cerr << "Unable to open file!" << endl;
        return 1;
    }

    for (int i = 0; i < num_of_maturity; ++i) {
        for (int j = 0; j < num_of_K; ++j) {
            file_sabr << out_sabr[0][j][i];
            if (j < num_of_K) file_sabr << ",";
        }
        file_sabr << "\n";
        for (int j = 0; j < num_of_K; ++j) {
            file_sabr << out_sabr[1][j][i];
            if (j < num_of_K) file_sabr << ",";
        }
        file_sabr << "\n";
    }
    file_sabr.close();
    std::cout << "Simulation result for SABR has been saved as CSV file!" << endl;
    std::cout << endl;

    // Memory Clear
    for (int i = 0; i < num_of_model; ++i) {
        for (int j = 0; j < num_of_K; ++j) {
            delete[] out_sabr[i][j];
        }
        delete[] out_sabr[i];
    }
    delete[] out_sabr;
    // First Experiment is terminated.


    return 0 ;
}

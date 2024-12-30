// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>

using namespace Rcpp;
using namespace Eigen;

class Cmpp {
private:
    Eigen::MatrixXd features;  // Matrix to store feature data
    Eigen::VectorXd x;         // Vector to store failure times
    Eigen::VectorXi delta1, delta2; // Vectors to store censoring indicators
    int nsamp; 
    int nfeatures; 
    double h;

public:
    // Constructor to initialize the class with data
    Cmpp(const Eigen::MatrixXd& features_, const Eigen::VectorXd& x_, const Eigen::VectorXi& delta1_, const Eigen::VectorXi& delta2_, const double h_) {
        features = features_; // Store the feature data as an Eigen matrix
        x = x_;               // Store the failure times as an Eigen vector
        delta1 = delta1_;     // Store censoring indicator 1 as an Eigen vector
        delta2 = delta2_;     // Store censoring indicator 2 as an Eigen vector
        nsamp = features.rows(); // number of samples 
        nfeatures = features.cols(); // number of features
        h = h_;
    }

    // Method to return features matrix
    Eigen::MatrixXd get_features() {
        return features;
    }

    // Method to return failure times (x)
    Eigen::VectorXd get_failure_times() {
        return x;
    }

     F_Gomp(double x, double alpha, double beta) {
        return 1 - std::exp(beta * (1 - std::exp(alpha * x)) / alpha); 
    }

     double f_Gomp(double x, double alpha, double beta) {
        return beta * std::exp(alpha * x + (beta / alpha) * (1 - std::exp(alpha * x)));
    }

    Rcpp::List GetNum_method() {
        return Rcpp::List::create(Named("Nsamp") = nsamp, 
                                  Named("Nfeature") = nfeatures);
    }

    double LogLike1_method(const Eigen::VectorXd& Param) {
        double likelihood = 0.0;
        const double epsilon = 1e-3; // Minimum threshold for stability

        for (int i = 0; i < nsamp; ++i) {
            double F1 = F_Gomp(x[i], Param[0], Param[1]);
            double F2 = F_Gomp(x[i], Param[2], Param[3]);
            double f1 = f_Gomp(x[i], Param[0], Param[1]);
            double f2 = f_Gomp(x[i], Param[2], Param[3]);

            double S = 1 - F1 - F2; // Overall survival function
            // Safeguard to handle numerical issues
            if (S <= 0) {
                S = epsilon;
            }

            // Compute log-likelihood based on delta indicators
            if (delta1[i] == 1) {
                likelihood += std::log(f1);
            } else if (delta2[i] == 1) {
                likelihood += std::log(f2);
            } else {
                likelihood += std::log(S);
            }
        }

        return -likelihood;
    }

    Eigen::VectorXd compute_numeric_gradient(const Eigen::VectorXd& Param) {
        Eigen::VectorXd grad = Eigen::VectorXd::Zero(Param.size()); 

        for (int i = 0; i < Param.size(); ++i) {
            Eigen::VectorXd Param_plus = Param;
            Eigen::VectorXd Param_minus = Param;

            Param_plus[i] += h; 
            Param_minus[i] -= h;

            double likelihood_plus = LogLike1_method(Param_plus);
            double likelihood_minus = LogLike1_method(Param_minus);

            grad[i] = (likelihood_plus - likelihood_minus) / (2 * h); 
        }

        return grad;
    }

    // Evaluate LogLik and Gradient
    double evaluate(const Eigen::VectorXd& Param, Eigen::VectorXd& grad) {
        grad = compute_numeric_gradient(Param);
        return -LogLike1_method(Param);
    }

    // Method to return censoring indicators (delta1 and delta2)
    Rcpp::List get_censoring_indicators() {
        return Rcpp::List::create(Named("delta1") = delta1, Named("delta2") = delta2);
    }
};

Cmpp* cmpp = nullptr; // Define a pointer to Cmpp class instance

// [[Rcpp::export]]
void Initialize(NumericMatrix features, NumericVector x, IntegerVector delta1, IntegerVector delta2, double h) {
    // Convert Rcpp types to Eigen types
    Eigen::Map<Eigen::MatrixXd> feature_matrix(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(features));  // Convert R matrix to Eigen matrix
    Eigen::Map<Eigen::VectorXd> x_vector(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(x));                // Convert R vector to Eigen vector
    Eigen::Map<Eigen::VectorXi> delta1_vector(Rcpp::as<Eigen::Map<Eigen::VectorXi>>(delta1));      // Convert R vector to Eigen vector
    Eigen::Map<Eigen::VectorXi> delta2_vector(Rcpp::as<Eigen::Map<Eigen::VectorXi>>(delta2));      // Convert R vector to Eigen vector

    // If an instance already exists, delete it
    if (cmpp != nullptr) {
        delete cmpp;
    }
    // Create a new instance of the Cmpp class
    cmpp = new Cmpp(feature_matrix, x_vector, delta1_vector, delta2_vector, h);
}

// [[Rcpp::export]]
double cdf_gomp(double time, double shape, double scale) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    return cmpp->F_Gomp(time, shape, scale); 
}

// [[Rcpp::export]]
double pdf_gomp(double x, double alpha, double beta) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    return cmpp->f_Gomp(x, alpha, beta); 
}

// [[Rcpp::export]]
Rcpp::List GetDim() {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    return cmpp->GetNum_method();
}

// [[Rcpp::export]]
SEXP LogLike1(SEXP paramSEXP) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }

    Eigen::Map<Eigen::VectorXd> Param(as<Eigen::Map<Eigen::VectorXd>>(paramSEXP));
    double result = cmpp->LogLike1_method(Param);
    return Rcpp::wrap(result);
}



// [[Rcpp::export]]
SEXP compute_grad(SEXP paramSEXP) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }

    Eigen::Map<Eigen::VectorXd> Param(as<Eigen::Map<Eigen::VectorXd>>(paramSEXP));
    Eigen::VectorXd grad = cmpp->compute_numeric_gradient(Param);
    return Rcpp::wrap(grad);
}


// [[Rcpp::export]]
Eigen::MatrixXd makeMat(int n, int m, double value){
    Eigen::MatrixXd mat = Eigen::MatrixXd::Constant(n, m, value);
    return mat; 
}

// Clean up memory by deleting the pointer when done
// [[Rcpp::export]]
void Cleanup() {
    if (cmpp != nullptr) {
        delete cmpp;
        cmpp = nullptr;
    }
}

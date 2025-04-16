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

    // for check length of vectors 
    void check_params_length(const Eigen::VectorXd& Params, int expected_length) {
        if (Params.size() != expected_length) {
            Rcpp::stop("The length of Params is incorrect. Expected length: " + std::to_string(expected_length) + ", but got: " + std::to_string(Params.size()));
        }
    }

    // Method to return features matrix
    Eigen::MatrixXd get_features() {
        return features;
    }

    // Method to return failure times (x)
    Eigen::VectorXd get_failure_times() {
        return x;
    }

    double F_Gomp(double x, double alpha, double beta) {
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
        int Len = 4; 
        check_params_length(Param, Len);
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

    // Method to compute the Hessian matrix
    Eigen::MatrixXd compute_numeric_hessian(const Eigen::VectorXd& Param) {
        int n = Param.size();
        Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(n, n);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Eigen::VectorXd Param_ij = Param;
                Param_ij[i] += h;
                Param_ij[j] += h;
                double f_ij = LogLike1_method(Param_ij);

                Param_ij[j] -= 2 * h;
                double f_i_j = LogLike1_method(Param_ij);

                Param_ij[i] -= 2 * h;
                double f_ij_ = LogLike1_method(Param_ij);

                Param_ij[j] += 2 * h;
                double f_i_j_ = LogLike1_method(Param_ij);

                hessian(i, j) = (f_ij - f_i_j - f_ij_ + f_i_j_) / (4 * h * h);
            }
        }

        return hessian;
    }

    // Method to compute F_cdf
    double F_cdf(const Eigen::VectorXd& Params, const Eigen::VectorXd& Z, double x) {
        double alpha = Params[0];
        double tau = Params[1];
        double rho = Params[2];
        rho = (rho < 0) * rho - (rho > 0) * rho;
        Eigen::VectorXd Beta = Params.tail(Params.size() - 3);
        double tempval = Z.dot(Beta);
        double temp1 = 1 - std::pow(1 + alpha * std::exp(tempval) * tau * (std::exp(rho * x) - 1) / rho, -1 / alpha);
        return temp1;
    }
    // Method to compute f_pdf
    double f_pdf(const Eigen::VectorXd& Params, const Eigen::VectorXd& Z, double x) {
        double alpha = Params[0];
        double tau = Params[1];
        double rho = Params[2];
        rho = (rho < 0) * rho - (rho > 0) * rho;
        Eigen::VectorXd Beta = Params.tail(Params.size() - 3);
        double tempval = Z.dot(Beta);
        double result = (tau * std::pow(alpha * tau * (std::exp(rho * x) - 1) * std::exp(tempval) / rho + 1, -1 / alpha) * 
                        std::exp(tempval) * std::exp(rho * x)) / 
                        (alpha * tau * (std::exp(rho * x) - 1) * std::exp(tempval) / rho + 1);
        return result;
    }
    double log_f(const Eigen::VectorXd& Params, const Eigen::MatrixXd& covars, const Eigen::VectorXd& x, const Eigen::MatrixXi& dData, int nk) {
        int n = covars.rows();
        int P = covars.cols();
        int ind = 3 + P;
        double s = 0.0;
        double Ftemp = 0.0;
        double fTemp = 0.0;
        for (int j = 0; j < n; ++j) {
            Ftemp = 0.0;
            fTemp = 0.0;
            for (int k = 0; k < nk; ++k) {
                Eigen::VectorXd parr = Params.segment(k * ind, ind);
                int delta = dData(j, k);
                fTemp += delta * std::log(f_pdf(parr, covars.row(j), x[j]));
                Ftemp += F_cdf(parr, covars.row(j), x[j]);
            }
            double Delta_Temp = 1.0 - dData.row(j).sum();
            if (Ftemp >= 1.0) {
                Ftemp = 0.5;
            }
            double F_final_Temp = Delta_Temp * std::log(1.0 - Ftemp);
            s += fTemp + F_final_Temp;
        }
        double tempres1 = s * s;
        return tempres1;
    }

    double log_f_single(const Eigen::VectorXd& Params) {
        int Len = 2 * (nfeatures + 3);
        check_params_length(Params, Len);
        Eigen::MatrixXi dData(nsamp, 2);
        dData.col(0) = delta1;
        dData.col(1) = delta2;
        return log_f(Params, features, x, dData, 2);
    }

    Eigen::VectorXd compute_log_f_gradient(const Eigen::VectorXd& Params) {
        Eigen::VectorXd grad = Eigen::VectorXd::Zero(Params.size());

        for (int i = 0; i < Params.size(); ++i) {
            Eigen::VectorXd Params_plus = Params;
            Eigen::VectorXd Params_minus = Params;

            Params_plus[i] += h;
            Params_minus[i] -= h;

            double log_f_plus = log_f_single(Params_plus);
            double log_f_minus = log_f_single(Params_minus);

            grad[i] = (log_f_plus - log_f_minus) / (2 * h);
        }

        return grad;
    }

    Eigen::MatrixXd compute_log_f_hessian(const Eigen::VectorXd& Params) {
        int n = Params.size();
        Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(n, n);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Eigen::VectorXd Params_ij = Params;
                Params_ij[i] += h;
                Params_ij[j] += h;
                double f_ij = log_f_single(Params_ij);

                Params_ij[j] -= 2 * h;
                double f_i_j = log_f_single(Params_ij);

                Params_ij[i] -= 2 * h;
                double f_ij_ = log_f_single(Params_ij);

                Params_ij[j] += 2 * h;
                double f_i_j_ = log_f_single(Params_ij);

                hessian(i, j) = (f_ij - f_i_j - f_ij_ + f_i_j_) / (4 * h * h);
            }
        }

        return hessian;
    }

// add 2025-03-27 
    double F_cdf2(const Eigen::VectorXd& Params, const Eigen::VectorXd& Z, double x){
        double tau = Params[0];
        double rho = Params[1];
        rho = (rho < 0) * rho - (rho > 0) * rho;
        Eigen::VectorXd Beta = Params.tail(Params.size() - 2);
        double tempval = Z.dot(Beta);
        double temp1 = tau * (std::exp(rho * x) - 1) * std::exp(tempval) / (tau * (std::exp(rho * x) - 1) * std::exp(tempval) + rho);
        return temp1;
    }

    double F_cdf3(const Eigen::VectorXd& Params, const Eigen::VectorXd& Z, double x){
        double tau = Params[0];
        double rho = Params[1];
        rho = (rho < 0) * rho - (rho > 0) * rho;
        Eigen::VectorXd Beta = Params.tail(Params.size() - 2);
        double tempval = Z.dot(Beta);
        double temp1 = 1 - std::exp(-(tau * (std::exp(rho * x) - 1) * std::exp(tempval))/rho);
        return temp1;
    }

    double f_pdf2(const Eigen::VectorXd& Params, const Eigen::VectorXd& Z, double x) {
        double tau = Params[0];
        double rho = Params[1];
        rho = (rho < 0) * rho - (rho > 0) * rho;
        Eigen::VectorXd Beta = Params.tail(Params.size() - 2);
        double tempval = Z.dot(Beta);
        double result1n = -std::pow(tau, 2)*(std::exp(rho * x)-1)*std::exp(rho*x)*std::exp(2 * tempval);
        double result1d = rho * std::pow(tau *(std::exp(rho*x)-1)*std::exp(tempval)/rho + 1, 2);
        double result2n = tau * std::exp(rho * x) * std::exp(tempval);
        double result2d = tau * (std::exp(rho * x) - 1)*std::exp(tempval)/rho + 1;
        return result1n/result1d + result2n/result2d;
    }
    double f_pdf3(const Eigen::VectorXd& Params, const Eigen::VectorXd& Z, double x) {
        double tau = Params[0];
        double rho = Params[1];
        rho = (rho < 0) * rho - (rho > 0) * rho;
        Eigen::VectorXd Beta = Params.tail(Params.size() - 2);
        double tempval = Z.dot(Beta);
        double result = tau * std::exp(rho * x) * (std::exp(-tau *(std::exp(rho * x) - 1)*std::exp(tempval)/rho)) * std::exp(tempval);
        return result;
    }
// end add 2025-03-27

// add in 2025-03-28
  double log_f2(const Eigen::VectorXd& Params, const Eigen::MatrixXd& covars, const Eigen::VectorXd& x, const Eigen::MatrixXi& dData, int nk) {
        int n = covars.rows();
        int P = covars.cols();
        int ind = 2 + P;
        double s = 0.0;
        double Ftemp = 0.0;
        double fTemp = 0.0;
        for (int j = 0; j < n; ++j) {
            Ftemp = 0.0;
            fTemp = 0.0;
            for (int k = 0; k < nk; ++k) {
                Eigen::VectorXd parr = Params.segment(k * ind, ind);
                int delta = dData(j, k);
                fTemp += delta * std::log(f_pdf2(parr, covars.row(j), x[j]));
                Ftemp += F_cdf2(parr, covars.row(j), x[j]);
            }
            double Delta_Temp = 1.0 - dData.row(j).sum();
            if (Ftemp >= 1.0) {
                Ftemp = 0.5;
            }
            double F_final_Temp = Delta_Temp * std::log(1.0 - Ftemp);
            s += fTemp + F_final_Temp;
        }
        double tempres1 = s * s;
        return tempres1;
    }

    double log_f_single2(const Eigen::VectorXd& Params) {
        int Len = 2 * (nfeatures + 2);
        check_params_length(Params, Len);
        Eigen::MatrixXi dData(nsamp, 2);
        dData.col(0) = delta1;
        dData.col(1) = delta2;
        return log_f2(Params, features, x, dData, 2);
    }

    Eigen::VectorXd compute_log_f_gradient2(const Eigen::VectorXd& Params) {
        Eigen::VectorXd grad = Eigen::VectorXd::Zero(Params.size());

        for (int i = 0; i < Params.size(); ++i) {
            Eigen::VectorXd Params_plus = Params;
            Eigen::VectorXd Params_minus = Params;

            Params_plus[i] += h;
            Params_minus[i] -= h;

            double log_f_plus = log_f_single2(Params_plus);
            double log_f_minus = log_f_single2(Params_minus);

            grad[i] = (log_f_plus - log_f_minus) / (2 * h);
        }

        return grad;
    }



  double log_f3(const Eigen::VectorXd& Params, const Eigen::MatrixXd& covars, const Eigen::VectorXd& x, const Eigen::MatrixXi& dData, int nk) {
        int n = covars.rows();
        int P = covars.cols();
        int ind = 2 + P;
        double s = 0.0;
        double Ftemp = 0.0;
        double fTemp = 0.0;
        for (int j = 0; j < n; ++j) {
            Ftemp = 0.0;
            fTemp = 0.0;
            for (int k = 0; k < nk; ++k) {
                Eigen::VectorXd parr = Params.segment(k * ind, ind);
                int delta = dData(j, k);
                fTemp += delta * std::log(f_pdf3(parr, covars.row(j), x[j]));
                Ftemp += F_cdf3(parr, covars.row(j), x[j]);
            }
            double Delta_Temp = 1.0 - dData.row(j).sum();
            if (Ftemp >= 1.0) {
                Ftemp = 0.5;
            }
            double F_final_Temp = Delta_Temp * std::log(1.0 - Ftemp);
            s += fTemp + F_final_Temp;
        }
        double tempres1 = s * s;
        return tempres1;
    }

    double log_f_single3(const Eigen::VectorXd& Params) {
        int Len = 2 * (nfeatures + 2);
        check_params_length(Params, Len);
        Eigen::MatrixXi dData(nsamp, 2);
        dData.col(0) = delta1;
        dData.col(1) = delta2;
        return log_f3(Params, features, x, dData, 2);
    }

    Eigen::VectorXd compute_log_f_gradient3(const Eigen::VectorXd& Params) {
        Eigen::VectorXd grad = Eigen::VectorXd::Zero(Params.size());

        for (int i = 0; i < Params.size(); ++i) {
            Eigen::VectorXd Params_plus = Params;
            Eigen::VectorXd Params_minus = Params;

            Params_plus[i] += h;
            Params_minus[i] -= h;

            double log_f_plus = log_f_single3(Params_plus);
            double log_f_minus = log_f_single3(Params_minus);

            grad[i] = (log_f_plus - log_f_minus) / (2 * h);
        }

        return grad;
    }

// end add in 2025-03-28

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
double cdf_gomp(double x, double alpha, double beta) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    return cmpp->F_Gomp(x, alpha, beta); 
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
SEXP LogLike1(SEXP param) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }

    Eigen::Map<Eigen::VectorXd> Param(as<Eigen::Map<Eigen::VectorXd>>(param));
    double result = cmpp->LogLike1_method(Param);
    return Rcpp::wrap(result);
}

// [[Rcpp::export]]
SEXP compute_grad(SEXP param) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }

    Eigen::Map<Eigen::VectorXd> Param(as<Eigen::Map<Eigen::VectorXd>>(param));
    Eigen::VectorXd grad = cmpp->compute_numeric_gradient(Param);
    return Rcpp::wrap(grad);
}

// [[Rcpp::export]]
SEXP compute_hessian(SEXP param) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }

    Eigen::Map<Eigen::VectorXd> Param(as<Eigen::Map<Eigen::VectorXd>>(param));
    Eigen::MatrixXd hessian = cmpp->compute_numeric_hessian(Param);
    return Rcpp::wrap(hessian);
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

// [[Rcpp::export]]
List bootstrap_variance(NumericMatrix features, NumericVector x, IntegerVector delta1, IntegerVector delta2, NumericVector initial_params, int n_bootstrap, std::string optimMethod) {
  int n = features.nrow();  // Number of rows in features
  int p = initial_params.size();  // Number of parameters
  Eigen::MatrixXd bootstrap_estimates(n_bootstrap, p);

  for (int i = 0; i < n_bootstrap; ++i) {
    // Sample with replacement (adjust for 0-based indexing)
    IntegerVector indices = sample(n, n, true) - 1;  // Subtract 1 to convert to 0-based indexing
    NumericMatrix features_boot(n, features.ncol());
    NumericVector x_boot(n);
    IntegerVector delta1_boot(n);
    IntegerVector delta2_boot(n);

    // Create bootstrap samples
    for (int j = 0; j < n; ++j) {
      features_boot(j, _) = features(indices[j], _);  // Access rows using adjusted indices
      x_boot[j] = x[indices[j]];
      delta1_boot[j] = delta1[indices[j]];
      delta2_boot[j] = delta2[indices[j]];
    }

    // Initialize with bootstrap sample
    Initialize(features_boot, x_boot, delta1_boot, delta2_boot, 1e-5);

    // Estimate parameters using optim from stats package
    Environment stats = Environment::namespace_env("stats");
    Function optim = stats["optim"];
    List optim_result = optim(
      Named("par") = initial_params,
      Named("fn") = Rcpp::InternalFunction(&LogLike1),
      Named("gr") = Rcpp::InternalFunction(&compute_grad),
      Named("method") = optimMethod
    );

    // Store the results in the bootstrap estimates matrix
    Eigen::VectorXd par = as<Eigen::VectorXd>(optim_result["par"]);
    bootstrap_estimates.row(i) = par;
  }

  // Compute variance of the bootstrap estimates
  Eigen::VectorXd variances(p);
  for (int j = 0; j < p; ++j) {
    variances[j] = (bootstrap_estimates.col(j).array() - bootstrap_estimates.col(j).mean()).square().sum() / (n_bootstrap - 1);
  }

  return List::create(
    Named("variances") = variances,
    Named("bootstrap_estimates") = bootstrap_estimates
  );
}


// [[Rcpp::export]]
double F_cdf_rcpp(NumericVector Params, NumericVector Z, double x) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    Eigen::VectorXd Params_eigen = as<Eigen::VectorXd>(Params);
    Eigen::VectorXd Z_eigen = as<Eigen::VectorXd>(Z);
    return cmpp->F_cdf(Params_eigen, Z_eigen, x);
}

// [[Rcpp::export]]
double f_pdf_rcpp(NumericVector Params, NumericVector Z, double x) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    Eigen::VectorXd Params_eigen = as<Eigen::VectorXd>(Params);
    Eigen::VectorXd Z_eigen = as<Eigen::VectorXd>(Z);
    return cmpp->f_pdf(Params_eigen, Z_eigen, x);
}

// [[Rcpp::export]]
double log_f_rcpp(NumericVector Params) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    Eigen::VectorXd Params_eigen = as<Eigen::VectorXd>(Params);
    return cmpp->log_f_single(Params_eigen);
}

// [[Rcpp::export]]
NumericVector compute_log_f_gradient_rcpp(NumericVector Params) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    Eigen::VectorXd Params_eigen = as<Eigen::VectorXd>(Params);
    Eigen::VectorXd grad = cmpp->compute_log_f_gradient(Params_eigen);
    return wrap(grad);
}

// [[Rcpp::export]]
NumericMatrix compute_log_f_hessian_rcpp(NumericVector Params) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    Eigen::VectorXd Params_eigen = as<Eigen::VectorXd>(Params);
    Eigen::MatrixXd hessian = cmpp->compute_log_f_hessian(Params_eigen);
    return wrap(hessian);
}

// add in 2025-03-28


// [[Rcpp::export]]
double F_cdf_rcpp2(NumericVector Params, NumericVector Z, double x) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    Eigen::VectorXd Params_eigen = as<Eigen::VectorXd>(Params);
    Eigen::VectorXd Z_eigen = as<Eigen::VectorXd>(Z);
    return cmpp->F_cdf2(Params_eigen, Z_eigen, x);
}

// [[Rcpp::export]]
double f_pdf_rcpp2(NumericVector Params, NumericVector Z, double x) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    Eigen::VectorXd Params_eigen = as<Eigen::VectorXd>(Params);
    Eigen::VectorXd Z_eigen = as<Eigen::VectorXd>(Z);
    return cmpp->f_pdf2(Params_eigen, Z_eigen, x);
}

// [[Rcpp::export]]
double log_f_rcpp2(NumericVector Params) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    Eigen::VectorXd Params_eigen = as<Eigen::VectorXd>(Params);
    return cmpp->log_f_single2(Params_eigen);
}

// [[Rcpp::export]]
NumericVector compute_log_f_gradient_rcpp2(NumericVector Params) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    Eigen::VectorXd Params_eigen = as<Eigen::VectorXd>(Params);
    Eigen::VectorXd grad = cmpp->compute_log_f_gradient2(Params_eigen);
    return wrap(grad);
}

// [[Rcpp::export]]
double F_cdf_rcpp3(NumericVector Params, NumericVector Z, double x) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    Eigen::VectorXd Params_eigen = as<Eigen::VectorXd>(Params);
    Eigen::VectorXd Z_eigen = as<Eigen::VectorXd>(Z);
    return cmpp->F_cdf3(Params_eigen, Z_eigen, x);
}

// [[Rcpp::export]]
double f_pdf_rcpp3(NumericVector Params, NumericVector Z, double x) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    Eigen::VectorXd Params_eigen = as<Eigen::VectorXd>(Params);
    Eigen::VectorXd Z_eigen = as<Eigen::VectorXd>(Z);
    return cmpp->f_pdf3(Params_eigen, Z_eigen, x);
}

// [[Rcpp::export]]
double log_f_rcpp3(NumericVector Params) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    Eigen::VectorXd Params_eigen = as<Eigen::VectorXd>(Params);
    return cmpp->log_f_single3(Params_eigen);
}

// [[Rcpp::export]]
NumericVector compute_log_f_gradient_rcpp3(NumericVector Params) {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    Eigen::VectorXd Params_eigen = as<Eigen::VectorXd>(Params);
    Eigen::VectorXd grad = cmpp->compute_log_f_gradient3(Params_eigen);
    return wrap(grad);
}


// end in 2025-03-28


// Start in 2025-04-14
// [[Rcpp::export]]
Rcpp::List GetData() {
    if (cmpp == nullptr) {
        Rcpp::stop("The Cmpp object has not been initialized.");
    }
    return Rcpp::List::create(
        Rcpp::Named("features") = cmpp->get_features(),
        Rcpp::Named("time") = cmpp->get_failure_times(),
        Rcpp::Named("delta1") = cmpp->get_censoring_indicators()["delta1"],
        Rcpp::Named("delta2") = cmpp->get_censoring_indicators()["delta2"]
    );
}
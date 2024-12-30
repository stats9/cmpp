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

# EBestimator
Simulation study code and results for author's manusript *A Nonparametric Bayes Approach for Sparse Sequence Estimation*.

We compared the following 10 methods for Sparse Gaussian Sequence Model:
- Martin and Walker (2014)'s empirical Bayes estimator 
- Koenker and Mizera (2014)'s empirical Bayes estimator
- Johnstone and Silverman (2005)'s empirical Bayes estimator
- Hard thresholding estimator
- Soft thresholding estimator
- SURE estimator
- Abramovich et al. (2006)'s FDR based estimator with parameters q=0.01
- Abramovich et al. (2006)'s FDR based estimator with parameters q=0.1
- Abramovich et al. (2006)'s FDR based estimator with parameters q=4
- Proposed Dirichlet process mixture based nonparametric Bayes Estimator

We conducted four simulation studies. Details are summarized in authors' manuscript. In RData file MSE and MAE are summarized. We have set seed
for all the simulation studies.

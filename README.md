This code is for computing the MAP under the horseshoe+ prior for the high-dimensional proportional hazards model. The data simulation is done in file data_simulation_HS.R while all necessary functions are saved into source_code_HS.R. All computation is done in R. 

For generating different datasets, the seed, and other parameters of the data need to be changed in data_simulation_HS.R. The output includes RMSE for the proposed method and LASSO and feature selection metrics under the proposed method and LASSO. Note that we are using two techniques of feature selection after obtaining the MAP estimator. These methods are referred to as HS1 and HS2. 


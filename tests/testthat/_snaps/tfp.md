# TFP model

    Call:
    TFPmodel(tsl = tsList, cycle = "RAR2")
    
    	State space model object of class TFPmodel
    
    cycle 				RAR2
    trend 				DT
    cubs
      cycle lags 			0
      error term			iid normal
      exogenous variables		-
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	36
      period 			1985 - 2020
      frequency 			annual
    
    Object is a valid object of class TFPmodel.

# TFP MLE fit

    Call:
    fit.TFPmodel(model = model)
    
    	State space model object of class TFPmodel
    
    cycle 				RAR2
    trend 				DT
    cubs
      cycle lags 			0
      error term			iid normal
      exogenous variables		-
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	36
      period 			1985 - 2020
      frequency 			annual
    
    
    	Maximum likelihood estimation results
    
    cycle
             Coefficient Standard Error t-statistic  p-value
      cA           0.533          0.153        3.49 0.000485
      cSigma       0.862          0.226        3.82 0.000133
      cTau         8.709          3.270        2.66 0.007746
    
    trend
              Coefficient Standard Error t-statistic p-value
      tdOmega      0.5016         0.3853        1.30   0.193
      tdPhi        0.8944         0.1019        8.77   0.000
      tdSigma      0.0777         0.0668        1.16   0.245
    
    cubs
              Coefficient Standard Error t-statistic p-value
      cuC0           2.72          0.296        9.17 0.00000
      cuConst        0.21          0.915        0.23 0.81825
      cuSigma        1.53          0.561        2.73 0.00631
      RMSE: 2.89
      R2: 0.367
      Box-Ljung test: X-squared = 2.17, df = 7.2, p-value = 0.956
    
             loglik             AIC             BIC             HQC signal-to-noise 
         -120.28316       258.56631       272.81799       263.54053         0.09013 

# TFP bayesian fit

    Call:
    fit.TFPmodel(model = model, method = "bayesian", R = 1000, thin = 2, 
        MLEfit = f)
    
    	State space model object of class TFPmodel
    
    cycle 				RAR2
    trend 				DT
    cubs
      cycle lags 			0
      error term			iid normal
      exogenous variables		-
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	36
      period 			1985 - 2020
      frequency 			annual
    
    
    	MCMC estimation results
    
    cycle
              Mean Median 85% HPDI-LB 85% HPDI-UB
      cA     0.417  0.410       0.164       0.656
      cSigma 1.077  0.915       0.521       1.462
      cTau   8.079  7.589       3.086      12.394
    
    trend
                Mean Median 85% HPDI-LB 85% HPDI-UB
      tdOmega 0.0145 0.0143   -0.000777      0.0258
      tdPhi   0.9447 0.9628    0.900142      0.9898
      tdSigma 0.0519 0.0424    0.015013      0.0832
    
    cubs
                  Mean   Median 85% HPDI-LB 85% HPDI-UB
      cuC0     2.74864 2.712603      2.2410      3.1618
      cuConst -0.00115 0.000207     -0.0722      0.0693
      cuSigma  2.78277 2.611966      1.6156      3.7409
              MRMSE signal-to-noise 
            1.48578         0.04823 


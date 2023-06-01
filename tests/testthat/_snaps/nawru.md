# NAWRU model

    Call:
    NAWRUmodel(tsl = tsList, exoType = exoType)
    
    	State space model object of class NAWRUmodel
    
    cycle 				AR2
    trend 				RW2
    phillips curve
      type 				TKP
      cycle lags 			0
      error term			iid normal
      exogenous variables		pcddws
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	59
      period 			1962 - 2020
      frequency 			annual
    
    Object is a valid object of class NAWRUmodel.

# NAWRU MLE fit

    Call:
    fit.NAWRUmodel(model = model, parRestr = parRestr)
    
    	State space model object of class NAWRUmodel
    
    cycle 				AR2
    trend 				RW2
    phillips curve
      type 				TKP
      cycle lags 			0
      error term			iid normal
      exogenous variables		pcddws
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	59
      period 			1962 - 2020
      frequency 			annual
    
    
    	Maximum likelihood estimation results
    
    cycle
             Coefficient Standard Error t-statistic  p-value
      cPhi1        1.261         0.1269       9.941 0.00e+00
      cPhi2       -0.398         0.4729      -0.841 4.00e-01
      cSigma       0.210         0.0405       5.188 2.12e-07
    
    trend
              Coefficient Standard Error t-statistic p-value
      tdSigma      0.0019       0.000689        2.76 0.00578
    
    phillips curve
              Coefficient Standard Error t-statistic  p-value
      pcC0       -0.35277          0.149     -2.3742 1.76e-02
      pcConst     0.00524          0.206      0.0254 9.80e-01
      pcSigma     1.21631          0.227      5.3509 8.75e-08
      pcddws     98.46177          9.691     10.1606 0.00e+00
      RMSE: 1.16
      R2: 0.628
      Box-Ljung test: X-squared = 21, df = 10, p-value = 0.0213
    
             loglik             AIC             BIC             HQC signal-to-noise 
         -132.80070       281.60141       298.22171       288.08930         0.00905 

# NAWRU bayesian fit

    Call:
    fit.NAWRUmodel(model = model, method = "bayesian", R = 10000, 
        burnin = 3000, thin = 10, MLEfit = f)
    
    	State space model object of class NAWRUmodel
    
    cycle 				AR2
    trend 				RW2
    phillips curve
      type 				TKP
      cycle lags 			0
      error term			iid normal
      exogenous variables		pcddws
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	59
      period 			1962 - 2020
      frequency 			annual
    
    
    	MCMC estimation results
    
    cycle
               Mean Median 85% HPDI-LB 85% HPDI-UB
      cPhi1   0.645  0.645       0.493      0.7952
      cPhi2  -0.151 -0.149      -0.261     -0.0235
      cSigma  0.150  0.145       0.079      0.2056
    
    trend
                Mean Median 85% HPDI-LB 85% HPDI-UB
      tdSigma 0.0243 0.0175     0.00228      0.0428
    
    phillips curve
                 Mean  Median 85% HPDI-LB 85% HPDI-UB
      pcC0    -0.9398 -0.9051     -1.5232      -0.265
      pcConst -0.0196 -0.0156     -0.2160       0.183
      pcSigma  3.1074  3.0504      2.3565       3.802
      pcddws   4.8980  5.0559      0.0722       9.974
              MRMSE signal-to-noise 
              1.793           0.162 


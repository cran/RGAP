# Kuttner model

    Call:
    KuttnerModel(tsl = tsList)
    
    	State space model object of class KuttnerModel
    
    cycle 				AR2
    trend 				RW1
    inflation equation
      cycle lags 			1
      error term			ARMA(0,3)
      exogenous variables		gdpGL1
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	30
      period 			1991 - 2020
      frequency 			annual
    
    Object is a valid object of class KuttnerModel.

# Kuttner MLE fit

    Call:
    fit.KuttnerModel(model = model, parRestr = parRestr)
    
    	State space model object of class KuttnerModel
    
    cycle 				AR2
    trend 				RW1
    inflation equation
      cycle lags 			1
      error term			ARMA(0,3)
      exogenous variables		gdpGL1
    anchor
      value 			-
      horizon 			-
    dimensions
      number of observations	30
      period 			1991 - 2020
      frequency 			annual
    
    
    	Maximum likelihood estimation results
    
    cycle
             Coefficient Standard Error t-statistic p-value
      cPhi1        1.462           0.15       9.750 0.00000
      cPhi2       -0.534           1.43      -0.375 0.70799
      cSigma       1.803           0.65       2.771 0.00558
    
    trend
              Coefficient Standard Error t-statistic  p-value
      tSigma        0.354          0.318        1.11 2.65e-01
      tdConst       1.995          0.278        7.17 7.28e-13
    
    inflation equation
                    Coefficient Standard Error t-statistic  p-value
      inflC1           -0.00125        0.00814      -0.154 0.877910
      inflConst        -0.14771        0.09812      -1.505 0.132221
      inflErrGamma1    -0.22683        0.24143      -0.940 0.347462
      inflErrGamma2    -0.29116        0.14544      -2.002 0.045287
      inflErrGamma3    -0.48240        0.22927      -2.104 0.035370
      inflGdpGL1        0.07391        0.04591       1.610 0.107460
      inflSigma         0.37865        0.10319       3.669 0.000243
      RMSE: 0.637
      R2: 0.385
      Box-Ljung test: X-squared = 2.54, df = 6, p-value = 0.864
    
             loglik             AIC             BIC             HQC signal-to-noise 
           -82.5270        189.0541        205.8685        194.4331          0.1964 


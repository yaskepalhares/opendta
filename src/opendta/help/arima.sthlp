{smcl}
{title:Title}

{p 4 4 2}
{bf:arima} {hline 2} ARIMA, ARMAX and other dynamic regression models

{title:Syntax}

{p 8 16 2}
{cmd:arima} depvar [indepvars] [if] [in] [{cmd:,} {cmd:ar(}numlist{cmd:)} {cmd:ma(}numlist{cmd:)} {cmd:arima(}#p{cmd:,}#d{cmd:,}#q{cmd:)} {cmd:sarima(}#P{cmd:,}#D{cmd:,}#Q{cmd:,}#s{cmd:)} {cmd:mar(}numlist{cmd:,}#s{cmd:)} {cmd:mma(}numlist{cmd:,}#s{cmd:)} {cmd:condition} {cmd:diffuse} {cmd:noconstant} {cmd:vce(opg}|{cmd:oim}|{cmd:robust)} {cmd:level(}#{cmd:)} {cmd:nolog}]

{title:Description}

{pstd}
Fits a linear regression whose disturbance follows an ARMA process, by exact maximum likelihood with the Kalman filter (the initial state covariance is the stationary one). {cmd:arima(p,d,q)} differences the dependent and independent variables d times and includes AR lags 1 to p and MA lags 1 to q. The default variance is OPG, as in Stata. The data must be {cmd:tsset} and the sample may not contain gaps. Seasonal terms are given with {cmd:sarima(P,D,Q,s)}, {cmd:mar(numlist,s)} and {cmd:mma(numlist,s)}; {cmd:condition} uses conditional maximum likelihood and {cmd:diffuse} a diffuse initial state covariance (10^6, VERIFY).

{title:Postestimation}

{phang2}{cmd:predict} newvar [, {cmd:xb} {cmd:y} {cmd:residuals} {cmd:structural}]{p_end}

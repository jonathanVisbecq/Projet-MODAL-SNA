//------------------------------------------------------------------------------
// Calibration du pas de temps.
//------------------------------------------------------------------------------

// Paramètres
lambda = 0.5
mu = 0.51
tmax = 800
N = 100
S = 2000
nbSimulations = 400
alpha = 0.95
h = 0.06

t = linspace(0, tmax, tmax+1)
[E1, ValConf1] = esperanceEncombrement(lambda, mu, t, nbSimulations, alpha)
[E2, ValConf2] = esperanceDiscr(lambda, mu, t, h, nbSimulations, alpha)

T = ones(2,1)*t

plot2d(t, E1)
V1 = [E1-ValConf1; E1+ValConf1]
plot2d(T', V1')

plot2d(t, E2)
V2 = [E2-ValConf2; E2+ValConf2]
plot2d(T', V2')
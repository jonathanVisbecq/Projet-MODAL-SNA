//------------------------------------------------------------------------------
// Calibration du pas de temps.
//------------------------------------------------------------------------------

scf(2)

// Paramètres
lambda = 0.5
mu = 0.51
tmax = 200
N = 50
S = 2000
nbSimulations = 100
h = 0.05

// Comparaison de l'esperance de l'encombrement en fonction du temps pour les methodes
// de simulation avec et sans discretisation.
t = linspace(0, tmax, tmax+1)
[E1, ValConf1] = esperanceEncombrement(lambda, mu, t, nbSimulations)
[E2, ValConf2] = esperanceDiscr(lambda, mu, t, h, nbSimulations)

T = ones(2,1)*t

V1 = [E1-ValConf1; E1+ValConf1]
plot2d(T', V1', style=[color("red"), color("red")])

V2 = [E2-ValConf2; E2+ValConf2]
plot2d(T', V2', style=[color("blue"),color("blue")])

title('Evolution de l''esperance de l''encombrement mémoire en fonction du temps', 'fontname', 'Symbol')

legend("ee", "rr")


// Comparaison des esperances des temps de saturation pour les methodes de simulation 
// avec et sans discretisation







// Comparaison des probabilites de depassement memoire pour les methodes de simulation 
// avec et sans discretisation


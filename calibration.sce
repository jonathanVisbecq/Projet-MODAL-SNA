//------------------------------------------------------------------------------
// Calibration du pas de temps.
//------------------------------------------------------------------------------

stacksize(150000000)

// Paramètres
lambda = 0.52
mu = 0.5
tmax = 2000
N = 50
S = 400
nbSimulations = 100
h = 0.05

// Comparaison de l'esperance de l'encombrement en fonction du temps pour les methodes de simulation avec et sans discretisation.
t = linspace(0, tmax, tmax+1)
[E1, ValConf1] = esperanceEncombrement(lambda, mu, t, nbSimulations)
[E2, ValConf2] = esperanceDiscr(lambda, mu, t, h, nbSimulations)

T = ones(2,1)*t

V1 = [E1-ValConf1; E1+ValConf1]
plot2d(T', V1', style=[color("red"), color("red")])

V2 = [E2-ValConf2; E2+ValConf2]
plot2d(T', V2', style=[color("blue"),color("blue")])

title('Evolution de l''esperance de l''encombrement mémoire en fonction du temps')

legend("Sans discretisation", "", "avec discretisation", "")


// Comparaison des esperances des temps de saturation pour les methodes de simulation 
// avec et sans discretisation

//[ETn1, ValConf1] = espTpsSat(lambda, mu, nbSimulations, N)
//[ETn2, ValConf2] = espTpsSatDiscr(lambda, mu, h, nbSimulations, N)
//
//disp('----- Intervalles de confiance pour l''esperance du temps de saturation -----')
//disp('Sans discretisation')
//disp(ETn1-ValConf1, ETn1+ValConf1)
//disp('Avec discretisation')
//disp(ETn2-ValConf1, ETn2+ValConf2)

//
//
//// Comparaison des probabilites de depassement memoire pour les methodes de simulation 
//// avec et sans discretisation

//ps1 = probabiliteSat(lambda, mu, nbSimulations, N, S)
//ps2 = probabiliteSatDiscr(lambda, mu, h, nbSimulations, N, S)
//
//disp('----- Probabilite de depassement memoire -----')
//disp('Sans discretisation')
//disp(ps1)
//disp('Avec discretisation')
//disp(ps2)
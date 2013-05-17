//------------------------------------------------------------------------------
// D�termination d'un crit�re de stabilit� du processus par m�thode graphique
//------------------------------------------------------------------------------

// Param�tres TEST
lambdaMin = 0.3
lambdaMax = 0.7
step = 0.05
mu = 0.5
tmax = 2000
nbSimulations = 100

// Calcul de l'esp�rance de l'encombrement en fonction du temps
t = linspace(0, tmax, tmax/5) 
E = []

for lambda=lambdaMin:step:lambdaMax
    E = [E; esperanceEncombrement(lambda, mu, t, nbSimulations)]
end

// Affichage
T = ones(length(lambda), 1)*t
plot2d(T', E')

str = "lambda = "
leg = []
for i=lambdaMin:step:lambdaMax 
    leg = [leg, str + string(i)]
end

legend(leg)
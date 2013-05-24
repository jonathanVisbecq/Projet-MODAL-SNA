//------------------------------------------------------------------------------
// Détermination d'un critère de stabilité du processus par méthode graphique
//------------------------------------------------------------------------------

// Paramètres
lambdaMin = 0.4
lambdaMax = 0.5
step = 0.02
mu = 0.5
tmax = 4000
nbSimulations = 200

// Calcul de l'espérance de l'encombrement en fonction du temps
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

// Coefficiens de regression lineaire avec le temps de l'esperance de l'encombrement
[a, b, sig] = reglin(T(1,:), E)
disp(a)
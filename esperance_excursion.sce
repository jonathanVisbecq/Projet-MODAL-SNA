//------------------------------------------------------------------------------
// Espérance d'une excursion
//------------------------------------------------------------------------------

stacksize(150000000)

// Paramètres
lambda=0.45;
mu=0.55;

// Nombre de simulations a effectuer
nbSimulations = 1000;
// Valeur du pas de discrétisation
h = 0.05;
// Nombre de variables aléatoires à simuler à chaque fois que nécessaire
n = 10000;


// Estimation de l'espérance de la longueur de l'excursion
function [l, lConf]=longueurExcursion(lambda, mu, nbSimulations, h, n)
    Tps = []
    for i=1:nbSimulations
        X = 1;
        Tn = 0;
        while X>0
            i = 1;;
            T = h*grand(1, n, 'geom', (lambda+mu)*h)
            U =  grand(1, n, 'def');
            e = 1*(U<=lambda/(lambda+mu)) + (-1)*(U>lambda/(lambda+mu));
            while (X>0) & (i<=n)
                Tn = Tn + T(i);
                X = X + e(i);
                i = i+1;
            end
        end
        Tps = [Tps, Tn];
    end
    l = sum(Tps)/nbSimulations
    lConf = 1.96*sqrt(variance(Tps))/sqrt(nbSimulations)
endfunction

//Tps = longueurExcursion(lambda, mu, nbSimulations, h, n);
//disp('Estimation et intervalle de confiance à 95%');
//ValConf = 1.96*sqrt(variance(Tps, 'c'))/sqrt(nbSimulations);
//E = sum(Tps, 'c')/nbSimulations;
//disp(E+ValConf, E, E-ValConf);


// Calcul de la valeur exacte
function l=longueurExacte(lambda, mu, nbSimulations, h, n)
    k_max = 10000;
    p = lambda*mu;
    s = 0;
    for k=1:k_max
        s = s + p ;
        p = p * (lambda*mu) * (2*k-1) * (2*k) / (k^2);
    end
endfunction

//disp('Valeur exacte');
//disp(4*s/(lambda+mu));









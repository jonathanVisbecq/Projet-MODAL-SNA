//------------------------------------------------------------------------------
// Espérance d'une excursion
//------------------------------------------------------------------------------

stacksize(150000000)

// Paramètres
lambda=0.45;
mu=0.55;

// Nombre de simulations a effectuer
nbSimulations = 3000;
// Valeur du pas de discr�tisation

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
            //T = h*grand(1, n, 'geom', (lambda+mu)*h)
            T = grand(1, n, 'exp', 1/(lambda+mu))
            U =  grand(1, n, 'def');
            e = 1*(U<=lambda/(lambda+mu)) + (-1)*(U>lambda/(lambda+mu));
            while (X>0) & (i<=n)
                Tn = Tn + T(i);
                X = X + e(i);
                i = i+1;
            end
        end
        //Tps = [Tps, Tn+h];
        Tps = [Tps, Tn]
    end
    l = sum(Tps)/nbSimulations
    lConf = 1.96*sqrt(variance(Tps))/sqrt(nbSimulations)
endfunction

[l, lConf] = longueurExcursion(lambda, mu, nbSimulations, h, n);
disp('Estimation et intervalle de confiance � 95%');
disp(l+lConf, l, l-lConf);



// Calcul de la valeur exacte

function l=longueurExacte(lambda, mu, nbSimulations, h, n)
    k_max = 10000;
    p = lambda*mu;
    s = 0;
    for k=1:k_max
        s = s + p ;
        p = p * (lambda*mu) * (2*k-1) * (2*k) / (k^2);
    end
    l = 4*s/(lambda+mu)
endfunction

disp('Valeur exacte');
disp(longueurExacte(lambda, mu, nbSimulations, h, n));









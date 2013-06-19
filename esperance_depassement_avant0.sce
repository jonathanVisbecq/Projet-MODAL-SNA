//------------------------------------------------------------------------------
// Calcul de la probabilitÃƒÂ© de dÃƒÂ©passement de la mÃƒÂ©moire avant le premier passage
// et de l'espÃƒÂ©rance de la temps de dÃƒÂ©passement sachant qu'il se produit avant le 
// premier retour en 0. -- sans changement de probabilitÃƒÂ© --
//------------------------------------------------------------------------------


stacksize(150000000)


// Paramètres
lambda=0.40
mu=0.6
N = 10

// Nombre de simulations a effectuer
nbSimulations = 500
// Valeur du pas de discrétisation
h = 0.05
// Nombre de variables aléatoires à simuler à chaque fois que nécessaire
n = 100


function [p, e, eConf]=probaEspDepassement(lambda, mu, N, nbSimulations, h, n)
    Tps = []
    nb = 0
    for i=1:nbSimulations
        X = 1
        Tn = 0
        while (X>0) & (X<N)
            i = 1
            //T = h*grand(1, n, 'geom', (lambda+mu)*h)
            T = grand(1, n, 'exp', 1/(lambda+mu))
            U =  grand(1, n, 'def')
            e = 1*(U<=lambda/(lambda+mu)) + (-1)*(U>lambda/(lambda+mu))
            while (X>0) & (X<N) & (i<=n)
                Tn = Tn + T(i)
                X = X + e(i)
                i = i+1
            end
        end
        if X==N then
            nb = nb + 1
            //Tps = [Tps, Tn+h]
            Tps = [Tps, Tn]
        end
    end
    p = nb/nbSimulations
    e = mean(Tps)
    disp(length(Tps))
    eConf = 1.96*sqrt(variance(Tps))/sqrt(nbSimulations)
endfunction


// Test de la validité de la formule
//Ed = []
//Eu = []
//Fu = []
//Fd = []
//
//nbSimulations2 = 5*nbSimulations
//m = 15
//[l, lConf] = longueurExcursion(lambda, mu, nbSimulations2, h, n)
//for N=1:m
//    [p, e, eConf] = probaEspDepassement(lambda, mu, N, nbSimulations2, h, n)
//    //[ETn, ValConf]=espTpsSatDiscr(lambda, mu, h, nbSimulations, N, n)
//    [ETn, ValConf]=espTpsSat(lambda, mu, nbSimulations, N)
//    Eu = [Eu, ((1-p)/p)*(l+lConf) + 0]
//    Ed = [Ed, ((1-p)/p)*(l-lConf) + 0]
//    Fu = [Fu, ETn+ValConf]
//    Fd = [Fd, ETn-ValConf]
//end
//
//plot2d(1:m, Eu, style=5)
//plot2d(1:m, Ed, style=5)
//plot2d(1:m, Fu, style=4)
//plot2d(1:m, Fd, style=4)



//disp(((1-p)/p)*(l+lConf) + e + eConf, ((1-p)/p)*l + e, ((1-p)/p)*(l-lConf) + e - eConf)
//disp(ETn+ValConf, ETn, ETn-ValConf)

//disp('Nombre de rÃƒÂ©alisations de l''ÃƒÂ©vÃƒÂ¨nement:')
//disp(nb)
//
//disp('ProbabilitÃƒÂ© de dÃƒÂ©passement de la mÃƒÂ©moire avant retour en zÃƒÂ©ro:')
//disp(p)

//[p, e, eConf]=probaEspDepassement(lambda, mu, N, nbSimulations, h, n)
//disp("Espérance du temps de dépassement sachant qu''il à lieu pendant")
//disp("la première excursion:")
//disp(e+eConf, e, e-eConf)


//
//disp(p)
//disp(((mu/lambda)-1)/(((mu/lambda)^N)-1))
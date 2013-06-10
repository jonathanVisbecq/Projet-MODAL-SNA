//------------------------------------------------------------------------------
// Calcul de la probabilit� de d�passement de la m�moire avant le premier passage
// et de l'esp�rance de la temps de d�passement sachant qu'il se produit avant le 
// premier retour en 0. -- sans changement de probabilit� --
//------------------------------------------------------------------------------


stacksize(150000000)

// Param�tres
lambda=0.45
mu=0.55
N = 20

// Nombre de simulations a effectuer
nbSimulations = 2000
// Valeur du pas de discr�tisation
h = 0.05
// Nombre de variables al�atoires � simuler � chaque fois que n�cessaire
n = 1000

function [p, e, eConf]=probaEspDepassement(lambda, mu, N, nbSimulations, h, n)
    Tps = []
    nb = 0
    for i=1:nbSimulations
        X = 1
        Tn = 0
        while (X>0) & (X<N)
            i = 1
            T = h*grand(1, n, 'geom', (lambda+mu)*h)
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
            Tps = [Tps, Tn]
        end
    end
    p = nb/nbSimulations
    e = sum(Tps)/nbSimulations
    eConf = 1.96*sqrt(variance(Tps))/sqrt(nbSimulations)
endfunction

[p, e, eConf] = probaEspDepassement(lambda, mu, N, nbSimulations, h, n)

//disp('Nombre de r�alisations de l''�v�nement:')
//disp(nb)
//
//disp('Probabilit� de d�passement de la m�moire avant retour en z�ro:')
//disp(p)
//
//disp("Esp�rance du temps de d�passement sachant qu''il � lieu pendant")
//disp("la premi�re excursion:")
//E = sum(Tps)/length(Tps)
//Valconf = 1.96 * sqrt(variance(Tps)) / sqrt(length(Tps))
//disp(E+Valconf, E, E-Valconf)

[l, lConf] = longueurExcursion(lambda, mu, nbSimulations, h, n)
p = ((mu/lambda)-1)/(((mu/lambda)^N)-1)
disp(((1-p)/p)*(l+lConf) + e + eConf, ((1-p)/p)*l + e, ((1-p)/p)*(l-lConf) + e - eConf)

[ETn, ValConf]=espTpsSatDiscr(lambda, mu, h, nbSimulations, N)
disp(ETn+ValConf, ETn, ETn-ValConf)
//
//disp(p)
//disp(((mu/lambda)-1)/(((mu/lambda)^N)-1))
//------------------------------------------------------------------------------
// Simule le syst�me avec discr�tisation en temps
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Simule une trajectoire de l'�tat du syst�me � l'aide d'une discr�tisation en
// temps.
//------------------------------------------------------------------------------
//
// lambda : (reel) Parametre de la loi d'arrivee des paquets
// mu : (reel) Parametre de la loi d'envoi des paquets.
// tmax : (reel) Temps jusqu'auquel on d�sire simuler le syst�me.
// h : (reel) Pas de la discr�tisation en temps.
// nbSimulations : (entier) Le nombre de simulations � effectuer    
//
// T : (vecteur ligne) Discr�tisation du temps entre 0 et tmax avec un pas h.
// X: (matrice) Valeur de l'encombrement aux instants donn�s par t. Chaque ligne 
//    correspond � une simulations.
//
function [T,X]=trajectoireDiscrete2(lambda, mu, tmax, h, nbSimulations)
    X = zeros(nbSimulations, 1)
    T = zeros(nbSimulations, 1)
    F = ones(nbSimulations, 1)
    Ones = ones(nbSimulations, 1)
    n = 1000
    while or(F)
        j = 1
        T1 = h*grand(nbSimulations, n, 'geom', lambda*h)
        T2 = h*grand(nbSimulations, n, 'geom', (lambda+mu)*h)
        U = grand(nbSimulations, n, 'def')
        e = 1*(U<=lambda/(lambda+mu)) + (-1)*(U>lambda/(lambda+mu))
        while or(F) & (j<=n)
            T = [T, T(:,$) + (T1(:,j).*(X(:,$)==0) + T2(:,j).*(X(:,$)>0))]
            X = [X, X(:,$) + ((Ones.*(X(:,$)==0)) + e(:,j).*(X(:,$)>0))]
            F = T(:,$)<tmax
            j = j+1
        end
    end
endfunction

[t,X] = trajectoireDiscrete2(0.36, 0.4,2000,1, 1)
plot2d(t,X, style=[color('red')])


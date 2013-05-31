//------------------------------------------------------------------------------
// Simule le syst�me avec discr�tisation en temps
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// Simule une trajectoire de l'�tat du syst�me � l'aide d'une discr�tisation en
// temps. ( /!\ Faux, ne plus utiliser /!\ )
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
function [T,X]=trajectoireDiscrete(lambda, mu, tmax, h, nbSimulations)
    imax = ceil(tmax/h)
    X = zeros(nbSimulations, 1)
    i = 1
    IncrA = 1*(grand(nbSimulations, imax+1, 'def')<=(lambda*h))
    IncrD = (-1)*(grand(nbSimulations, imax+1, 'def')<=(mu*h))
    while i<=imax
        X = [X, X(:,$) + IncrA(:,i).*(X(:,$)==0) + (IncrD(:,i) + IncrA(:,i)).*(X(:,$)>0)]
        i = i+1     
    end
    T = linspace(0, imax*h, imax+1)   
endfunction


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

[t,X] = trajectoireDiscrete2(0.45,0.5,1000,0.05, 1)
plot2d(t,X, style=[color('red')])


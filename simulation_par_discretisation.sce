//------------------------------------------------------------------------------
// Simule le syst�me avec discr�tisation en temps
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// Simule une trajectoire de l'�tat du syst�me � l'aide d'une discr�tisation en
// temps. (PREMIERE VERSION --> n'utilise pas de loi g�om�trique)
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

//[t,X] = trajectoireDiscrete(0.4,0.5,10000,1, 1)
//disp(size(t))
//disp(size(X))
//plot2d(t,X, style=[color('red')])


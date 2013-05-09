//------------------------------------------------------------------------------
// Simule le système avec discrétisation en temps
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// Simule une trajectoire de l'état du système à l'aide d'une discrétisation en
// temps. (PREMIERE VERSION --> n'utilise pas de loi géométrique)
//------------------------------------------------------------------------------
//
// lambda : (reel) Parametre de la loi d'arrivee des paquets
// mu : (reel) Parametre de la loi d'envoi des paquets.
// tmax : (reel) Temps jusqu'auquel on désire simuler le système.
// h : (reel) Pas de la discrétisation en temps.
//
// t : (vecteur ligne) Discrétisation du temps entre 0 et tmax avec un pas h.
// X: (vecteur ligne) Valeur de l'encombrement aux instants donnés par t.
//
function [t,X]=trajectoireDiscrete(lambda, mu, tmax, h)
    imax = ceil(tmax/h)
    X = [0]
    i = 2
    Ua = grand(1, imax+1, 'def')
    A = 1*(Ua<(lambda*h))
    Ud = grand(1, imax+1, 'def')
    D = (-1)*(Ud<=(mu*h))
    while i<=(imax+1)
        if(X(i-1)==0)
            X(i) = X(i-1) + A(i)
        else
            X(i) = X(i-1) + A(i) + D(i)
        end
        i = i+1     
    end
    t = linspace(0, imax*h, imax+1)   
endfunction

[t,X] = trajectoireDiscrete(0.3,0.5,1000,0.05)
plot2d(t,X, style=[color('red')])
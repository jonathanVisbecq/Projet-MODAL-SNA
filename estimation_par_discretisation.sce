//------------------------------------------------------------------------------
// Estimations utilisant la simulation par discr�tisation en temps
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Evolution de l'esp�rance de l'encombrement avec le temps
//------------------------------------------------------------------------------
//
// lambda (reel) : Parametre de la loi d'arrivee des paquets
// mu (reel) : Parametre de la loi de la duree d'envoie des paquets
// t (vecteur ligne) : Temps pour lesquels estimer l'esp�rance de l'encombrement.
// h : (reel) Pas de la discr�tisation en temps.
// nbsimulations : (entier) Le nombre de simulations utilis�es pour calculer
// l'esp�rance.
//
// E (vecteur ligne) : Les esp�rances correspondant aux temps t.
// ValConf (vecteur ligne) : L'ecart entre la moyenne empirique et les bornes
//                          de l'intervalle de confiance
//
function [E, ValConf]=esperanceDiscr(lambda, mu, t, h, nbSimulations)
    tmax = t($)
    [T,X] = trajectoireDiscrete(lambda, mu, tmax, h, nbSimulations)
    Y = []
    for i=1:nbSimulations
        [ind, occ, inf] = dsearch(t, T)
        Y = [Y; X(i, ind)]
    end
    ValConf = 1.96*sqrt(variance(Y, 'r'))/sqrt(nbSimulations)
    E = sum(Y, 'r')/nbSimulations
endfunction

//------------------------------------------------------------------------------
// Simulation du temps de saturation.
//------------------------------------------------------------------------------
//
// lambda (reel) : Parametre de la loi d'arrivee des paquets
// mu (reel) : Parametre de la loi de la duree d'envoie des paquets
// h : (reel) Pas de la discr�tisation en temps.
// nbSimulation (entier) : Nombre de simulations a effectuer.
// N (entier) : Taille de la m�moire tampon.
//
// Tn (vecteur ligne) : Les temps de saturation donn�es par les simulations
//
function Tn=tempsDeSaturationDiscr(lambda, mu, h, nbSimulations, N)
    n = 1000
    X = zeros(nbSimulations, 1)
    Tn = zeros(nbSimulations, 1)
    F = ones(nbSimulations, 1)
    Ones = ones(nbSimulations, 1)
    while or(F)
        i = 1
        IncrA = 1*(grand(nbSimulations, n, 'def')<=(lambda*h))
        IncrD = (-1)*(grand(nbSimulations, n, 'def')<=(mu*h))
        while or(F) & (i<=n)
            X = X + F.*(IncrA(:,i).*(X==0) + (IncrD(:,i) + IncrA(:,i)).*(X>0))         
            Tn = Tn + F.*Ones
            F = X<N
            i = i+1
        end
    end
    Tn = (Tn*h)'
endfunction

//------------------------------------------------------------------------------
// Estimation de la probabilit� de saturation avant un temps donn�.
//------------------------------------------------------------------------------
//
// lambda (reel) : Parametre de la loi d'arrivee des paquets
// mu (reel) : Parametre de la loi de la duree d'envoie des paquets
// h : (reel) Pas de la discr�tisation en temps.
// nbSimulation (entier) : Nombre de simulations a effectuer.
// N (entier) : Taille de la m�moire tampon.
// S (reel) : Le temps S dans la d�finition de la probabilit� (cf �nonc�)
//
// pS (reel) : Estimation de la probabilit� de saturation de la m�moire tampon
//             avant le temps S.
//
function pS=probabiliteSatDiscr(lambda, mu, h, nbSimulations, N, S)
    Tn = tempsDeSaturationDiscr(lambda, mu, h, nbSimulations, N) 
    pS = sum(Tn<=S)/nbSimulations   
endfunction


//n = 1000
//t = linspace(0, n, n+1)
//E = esperanceDiscr(0.5, 0.5, t, 1, 500)
//plot(t,E)
//
//Tn = tempsDeSaturationDiscr(0.55, 0.5, 0.1, 1000, 50)
//disp(sum(Tn)/1000)

//disp(probabiliteSatDiscr(0.6, 0.5, 0.1, 1000, 50, 400))
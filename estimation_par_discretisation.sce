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
//z
function [E, ValConf]=esperanceDiscr(lambda, mu, t, h, nbSimulations)
    tmax = t($)
    [T,X] = trajectoireDiscrete2(lambda, mu, tmax, h, nbSimulations)
    Y = []
    for i=1:nbSimulations
        [ind, occ, inf] = dsearch(t, T(i, :))
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
        T1 = h*grand(nbSimulations, n, 'geom', lambda*h)
        T2 = h*grand(nbSimulations, n, 'geom', (lambda+mu)*h)
        U = grand(nbSimulations, n, 'def')
        e = 1*(U<=lambda/(lambda+mu)) + (-1)*(U>lambda/(lambda+mu))
        while or(F) & (i<=n)
            Tn = Tn + F.*(T1(:,i).*(X==0) + T2(:,i).*(X>0))
            X = X + F.*((Ones.*(X==0)) + e(:,i).*(X>0))
            F = X<N
            i = i+1
        end
    end
    Tn = Tn'
endfunction


//------------------------------------------------------------------------------
// Estimation de l'esp�rance du temps de saturation de la memoire
//------------------------------------------------------------------------------
//
// lambda (reel) : Parametre de la loi d'arrivee des paquets
// mu (reel) : Parametre de la loi de la duree d'envoie des paquets
// h : (reel) Pas de la discr�tisation en temps.
// nbSimulation (entier) : Nombre de simulations a effectuer.
// N (entier) : Taille de la m�moire tampon.
//
// ETn (reel) : Approximation de l'esperance du temps de saturation
// ValConf (vecteur ligne) : L'ecart entre la moyenne empirique et les bornes
//                          de l'intervalle de confiance
//
function [ETn, ValConf]=espTpsSatDiscr(lambda, mu, h, nbSimulations, N)
    Tn = tempsDeSaturationDiscr(lambda, mu, h, nbSimulations, N)
    ETn = sum(Tn)/nbSimulations
    ValConf = 1.96*sqrt(variance(Tn, 'c'))/sqrt(nbSimulations)
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
//E = esperanceDiscr(0.45, 0.5, t, 0.01, 200)
//plot(t,E)
//

stacksize(150000000)
e = espTpsSatDiscr(0.2/0.05, 0.5/0.05, 0.05, 1000, 3)
disp(e)

//disp(probabiliteSatDiscr(0.6, 0.5, 0.1, 1000, 50, 400))n
//------------------------------------------------------------------------------
// Estimations basées sur les simulations exactes.
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// Evolution de l'espérance de l'encombrement avec le temps
//------------------------------------------------------------------------------
//
// lambda (reel) : Parametre de la loi d'arrivee des paquets
// mu (reel) : Parametre de la loi de la duree d'envoie des paquets
// t (vecteur ligne) : Temps pour lesquels estimer l'espérance de l'encombrement.
// nbsimulations : (entier) Le nombre de simulations utilisées pour calculer
// l'espérance.
//
// E (vecteur ligne) : Les espérances correspondant aux temps t.
//
function [E, ValConf]=esperanceEncombrement(lambda, mu, t, nbSimulations)
    X = trajectoireExacte1(lambda, mu, t, nbSimulations)
    ValConf = 1.96*sqrt(variance(X, 'r'))/sqrt(nbSimulations)
    E = sum(X, 'r')/nbSimulations
endfunction

//------------------------------------------------------------------------------
// Simulation du temps de saturation.
//------------------------------------------------------------------------------
//
// lambda (reel) : Parametre de la loi d'arrivee des paquets
// mu (reel) : Parametre de la loi de la duree d'envoie des paquets
// nbSimulation (entier) : Nombre de simulations a effectuer.
// N (entier) : Taille de la mémoire tampon.
//
// Tn (vecteur ligne) : Les temps de saturation données par les simulations
//
function Tn=tempsDeSaturation(lambda, mu, nbSimulations, N)
    n = 2000
    X = zeros(nbSimulations, 1)
    Tn = zeros(nbSimulations, 1)
    F = ones(nbSimulations, 1)
    while or(F)
        i = 1
        T1 = grand(nbSimulations, n, 'exp', 1/lambda)
        T2 = grand(nbSimulations, n, 'exp', 1/(lambda+mu))
        U = grand(nbSimulations, n, 'def')
        e = 1*(U<=lambda/(lambda+mu)) + (-1)*(U>lambda/(lambda+mu))
        while or(F) & (i<=n)
            Tn = Tn + F.*(T1(:,i).*(X==0) + T2(:,i).*(X>0))
            X = X + F.*(ones(nbSimulations, 1).*(X==0) + e(:,i).*(X>0))         
            F = X<N
            i = i+1
        end
    end
    Tn = Tn'
endfunction

//------------------------------------------------------------------------------
// Estimation de l'espérance du temps de saturation de la memoire
//------------------------------------------------------------------------------
//
// lambda (reel) : Parametre de la loi d'arrivee des paquets
// mu (reel) : Parametre de la loi de la duree d'envoie des paquets
// nbSimulation (entier) : Nombre de simulations a effectuer.
// N (entier) : Taille de la mémoire tampon.
//
// ETn (reel) : Approximation de l'esperance du temps de saturation
// ValConf (vecteur ligne) : L'ecart entre la moyenne empirique et les bornes
//                          de l'intervalle de confiance
//
function [ETn, ValConf]=espTpsSat(lambda, mu, nbSimulations, N)
    Tn=tempsDeSaturation(lambda, mu, nbSimulations, N)
    ETn = sum(Tn)/nbSimulations
    ValConf = 1.96*sqrt(variance(Tn, 'c'))/sqrt(nbSimulations)
endfunction



//------------------------------------------------------------------------------
// Estimation de la probabilité de saturation avant un temps donné.
//------------------------------------------------------------------------------
//
// lambda (reel) : Parametre de la loi d'arrivee des paquets
// mu (reel) : Parametre de la loi de la duree d'envoie des paquets
// nbSimulation (entier) : Nombre de simulations a effectuer.
// N (entier) : Taille de la mémoire tampon.
// S (reel) : Le temps S dans la définition de la probabilité (cf énoncé)
//
// pS (reel) : Estimation de la probabilité de saturation de la mémoire tampon
//             avant le temps S.
// ValConf (vecteur ligne) : L'ecart entre la moyenne empirique et les bornes
//                          de l'intervalle de confiance a 95%
//
function pS=probabiliteSat(lambda, mu, nbSimulations, N, S)
    Tn = tempsDeSaturation(lambda, mu, nbSimulations, N)
    pS = sum(Tn<=S)/nbSimulations
endfunction

//stacksize(150000000)
//n = 500
//t = linspace(0, n, n+1)
//E = esperanceEncombrement(0.45, 0.5, t, 1000)
//plot(t,E)


e = espTpsSat(0.45, 0.5, 1000, 3)
disp(e)


//disp(probabiliteSat(0.6, 0.5, 1000, 50, 400))









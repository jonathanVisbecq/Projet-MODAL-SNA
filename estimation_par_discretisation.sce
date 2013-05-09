//------------------------------------------------------------------------------
// Estimations utilisant la simulation par discrétisation en temps
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Evolution de l'espérance de l'encombrement avec le temps
//------------------------------------------------------------------------------
//
// lambda (reel) : Parametre de la loi d'arrivee des paquets
// mu (reel) : Parametre de la loi de la duree d'envoie des paquets
// t (vecteur ligne) : Temps pour lesquels estimer l'espérance de l'encombrement.
// h : (reel) Pas de la discrétisation en temps.
// nbsimulations : (entier) Le nombre de simulations utilisées pour calculer
// l'espérance.
//
// E (vecteur ligne) : Les espérances correspondant aux temps t.
//
function E=esperanceDiscr(lambda, mu, t, h, nbSimulations)
    tmax = t($)
    [T,X] = trajectoireDiscrete(lambda, mu, tmax, h, nbSimulations)
    Y = []
    for i=1:nbSimulations
        [ind, occ, inf] = dsearch(t, T)
        Y = [Y; X(i, ind)]
    end

    E = sum(Y, 'r')/nbSimulations
endfunction

//------------------------------------------------------------------------------
// Simulation du temps de saturation.
//------------------------------------------------------------------------------
//
// lambda (reel) : Parametre de la loi d'arrivee des paquets
// mu (reel) : Parametre de la loi de la duree d'envoie des paquets
// h : (reel) Pas de la discrétisation en temps.
// nbSimulation (entier) : Nombre de simulations a effectuer.
// N (entier) : Taille de la mémoire tampon.
//
// Tn (vecteur ligne) : Les temps de saturation données par les simulations
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
// Estimation de la probabilité de saturation avant un temps donné.
//------------------------------------------------------------------------------
//
// lambda (reel) : Parametre de la loi d'arrivee des paquets
// mu (reel) : Parametre de la loi de la duree d'envoie des paquets
// h : (reel) Pas de la discrétisation en temps.
// nbSimulation (entier) : Nombre de simulations a effectuer.
// N (entier) : Taille de la mémoire tampon.
// S (reel) : Le temps S dans la définition de la probabilité (cf énoncé)
//
// pS (reel) : Estimation de la probabilité de saturation de la mémoire tampon
//             avant le temps S.
//
function pS=probabiliteSatDiscr(lambda, mu, h, nbSimulations, N, S)
    Tn = tempsDeSaturationDiscr(lambda, mu, h, nbSimulations, N)
    pS = sum(Tn<S)/nbSimulations
endfunction






//n = 1000
//t = linspace(0, n, n+1)
//E = esperanceDiscr(0.5, 0.5, t, 1, 500)
//plot(t,E)
//
//Tn = tempsDeSaturationDiscr(0.55, 0.5, 0.1, 1000, 50)
//disp(sum(Tn)/1000)

disp(probabiliteSatDiscr(0.6, 0.5, 0.1, 1000, 50, 400))
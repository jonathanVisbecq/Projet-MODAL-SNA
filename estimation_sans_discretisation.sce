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
//
// E (vecteur ligne) : Les espérances correspondant aux temps t.
//
function E=esperanceEncombrement(lambda, mu, t)
    n = 1500
    X = trajectoireExacte1(lambda, mu, t, n)
    E = sum(X, 'r')/n
endfunction

//------------------------------------------------------------------------------
// Estimation du temps de saturation.
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
//
function pS=probabiliteSaturation(lambda, mu, nbSimulations, N, S)
    Tn = tempsDeSaturation(lambda, mu, nbSimulations, N)
    pS = sum(Tn<S)/nbSimulations
endfunction



//t = linspace(0, 1000, 1000)
//E = esperanceEncombrement(0.6, 0.5, t)
//plot(t,E)

//Tn = tempsDeSaturation(0.6, 0.5, 1000, 50)
//disp(sum(Tn)/1000)

//disp(probabiliteSaturation(0.6, 0.5, 5000, 50, 200))









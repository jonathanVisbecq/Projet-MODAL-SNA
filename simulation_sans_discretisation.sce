//------------------------------------------------------------------------------
// Simulation trajectoire exacte sans discr�tisation en temps.
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// Renvoie une matrice � deux lignes repr�sentant les temps d'arriv� et d'envoie
// des paquets
//------------------------------------------------------------------------------
//
// lambda : (reel) Parametre de la loi d'arrivee des paquets
// mu : (reel) Parametre de la loi d'envoi des paquets.
// tmax : (reel) Temps jusqu'auquel on veut pouvoir simuler le syst�me
// nbSimulations : (entier) Le nombre de simulations � effectuer
//
// T : (matrice) Le nombre de colonnes est le nombre de sauts n�cessaires pour
//    parvenir au temps tmax. Les lignes correspondent au diff�rentes simulations.
//    Les valeurs sont les temps successifs d'arriv� des paquets.`
// V : (matrice) Le nombre de colonnes est le nombre de sauts n�cessaires pour
//    parvenir au temps tmax. Les lignes correspondent au diff�rentes simulations.
//    Les valeurs sont les temps successifs de d�part des paquets.
//
function [T,V]=simulationExacte(lambda, mu, tmax, nbSimulations)
    T = zeros(nbSimulations, 1)
    T = []
    F = ones(nbSimulations, 1)
    n = 500
    // Calcul des temps d'arriv�e des paquets. On en ajoute tant que l'on n'a
    // pas atteint tmax.
    while or(F)
        Tinter = grand(nbSimulations, n, 'exp', 1/lambda)
        T = [T, cumsum(Tinter, 'c') + T(:,$)*ones(1,n)]
        F = T(:,$)<tmax
    end

    
    // Calcul des temps d'envoie de paquets. V_{i} est le temps de d�part du
    // paquet i
    E = grand(nbSimulations, length(T(1,:)), 'exp', 1/mu)
    V = [T(:,1)+E(:,1)]
    for i = 2:length(T(1,:))
        V = [V, max(V(:,i-1), T(:,i))+E(:,i)]
    end
endfunction

//------------------------------------------------------------------------------
// Simule une trajectoire de l'�tat du syst�me � l'aide de la fonction auxiliaire
// 'simulationExacte'
//------------------------------------------------------------------------------
//
// lambda : (reel) Parametre de la loi d'arrivee des paquets
// mu : (reel) Parametre de la loi d'envoi des paquets.
// x : (vecteur ligne) Temps pour lesquels on d�sire obtenir l'�tat du syst�me
// nbSimulations : (entier) Le nombre de simulations � effectuer
//
// X: (matrice) Valeur de l'encombrement aux instants donn�s par x. Chaque ligne
//   repr�sente une simulation.
//
function X=trajectoireExacte1(lambda, mu, x, nbSimulations)
    tmax = x($)
    X = []
    [T,V] = simulationExacte(lambda, mu, tmax, nbSimulations)
    for i=1:nbSimulations
        [indA, occA, infA] = dsearch(x, T(i,:))
        [indD, occD, infD] = dsearch(x, V(i,:))
        X = [X; indA - indD]
    end
endfunction


//------------------------------------------------------------------------------
// Simule une trajectoire de l'�tat du syst�me sans fonction auxiliaire par la 
// m�thode de l'�nonc�. Ici, on ne choisit pas les temps pour lesquels on calcule
// X: il s'agit des temps pour lesquels l'encombrement varie.
//------------------------------------------------------------------------------
//
// lambda : (reel) Parametre de la loi d'arrivee des paquets
// mu : (reel) Parametre de la loi d'envoi des paquets.
// tmax : (reel) Temps jusqu'auquel on veut pouvoir simuler le syst�me
// nbSimulations : (entier) Le nombre de simulations � effectuer
//
// T : (matrice) Temps pour lesquels la valeur de l'encombrement change. Chaque 
//    ligne repr�sente une simulation.
// X : (vecteur ligne) Valeurs correspondantes de l'encombrement. Chaque ligne
//    repr�sente une simulation.
//
function [T,X]=trajectoireExacte2(lambda, mu, tmax, nbSimulations)
    X = zeros(nbSimulations, 1)
    T = zeros(nbSimulations, 1)
    F = ones(nbSimulations, 1)
    Ones = ones(nbSimulations, 1)
    n = 1000
    while or(F)
        j = 1
        T1 = grand(nbSimulations, n, 'exp', 1/lambda)
        T2 = grand(nbSimulations, n, 'exp', 1/(lambda+mu))
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

function [T,X]=trajectoireExacte22(lambda, mu, n, nbSimulations)
    X = zeros(nbSimulations, 1)
    T = zeros(nbSimulations, 1)
    Ones = ones(nbSimulations, 1)
    T1 = grand(nbSimulations, n, 'exp', 1/lambda)
    T2 = grand(nbSimulations, n, 'exp', 1/(lambda+mu))
    U = grand(nbSimulations, n, 'def')
    e = 1*(U<=lambda/(lambda+mu)) + (-1)*(U>lambda/(lambda+mu))
    for j=1:n
        T = [T, T(:,$) + (T1(:,j).*(X(:,$)==0) + T2(:,j).*(X(:,$)>0))]
        X = [X, X(:,$) + ((Ones.*(X(:,$)==0)) + e(:,j).*(X(:,$)>0))]
    end
endfunction


//------------------------------------------------------------------------------
// Repr�sente une trajectoire de l'encombrement m�moire.
//------------------------------------------------------------------------------
//
// lambda : (reel) Parametre de la loi d'arrivee des paquets
// mu : (reel) Parametre de la loi d'envoi des paquets.
// tmax : (reel) Temps jusqu'auquel on veut pouvoir simuler le syst�me
// nbPoints: (vecteur ligne) Points en lesquels afficher l'�tat du syst�me, si 
//           l'on utilise 'trajectoireExacte1' pour la simulation
//
// traj (string) : (soit 'traj1' soit 'traj2') D�termine quel m�thode on applique pour
//        simuler la trajectoire
//
function evolutionTrajectoire(lambda, mu, tmax, nbPoints, traj)
    if traj=='traj1' then
        x = linspace(0, tmax, nbPoints)
        X = trajectoireExacte1(lambda, mu, x, 1)
    elseif traj=='traj2' then
        [x, X] = trajectoireExacte2(lambda, mu, tmax, 1)
    end
        
    plot2d2(x, X, style=[color('red')])
    legend('Evolution de l''encombrement m�moire en fonction du temps', pos=-5)
    l=get('current_entity')
    l.font_style = 5;
    l.font_size=2
    xlabel('Temps', 'fontsize', 2)
    ylabel('Encombrement', 'fontsize', 2)
endfunction













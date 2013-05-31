//------------------------------------------------------------------------------
// Changement de probabilité pour calculer E[Tn]
// 
//------------------------------------------------------------------------------

clear
stacksize(150000000)

N = 3
nbSimulations = 1000
h = 0.05

lambda=0.2/h
mu=0.5/h

new_lambda = mu
new_mu = lambda


function P=matriceTransition(a,b)
    P = zeros(N+1, N+1)
    P(N+1, N+1) = 1
    P(1,1) = 1-a*h
    P(1,2)= a*h
    for i=2:N
        P(i,i) = 1 - b*h - a*h
        P(i,i-1) = b*h
        P(i,i+1) = a*h
    end
endfunction



// Definition de P (matrice de transition du systeme discretise)
P = matriceTransition(lambda, mu)

// Definition de Q (matrice de transition du changement de probabilite)
Q = matriceTransition(new_lambda, new_mu)

//n = 1000
//X1 = ones(nbSimulations, 1)
//X2 = ones(nbSimulations, 1)
//Tn = zeros(nbSimulations, 1)
//F = ones(nbSimulations, 1)
//L = ones(nbSimulations, 1)
//Ones = ones(nbSimulations, 1)
//while or(F)
//    i = 1
//    T1 = grand(nbSimulations, n, 'geom', new_lambda*h)
//    T2 = grand(nbSimulations, n, 'geom', (new_lambda+new_mu)*h)
//    U = grand(nbSimulations, n, 'def')
//    e = 1*(U<=new_lambda/(new_lambda+new_mu)) + (-1)*(U>new_lambda/(new_lambda+new_mu))
//    while or(F) & (i<=n)
//        K = (T1(:,i).*(X1==1) + T2(:,i).*(X1>1))
//        Tn = Tn + h*F.*K
//        X2 = X1 + F.*((Ones.*(X1==1)) + e(:,i).*(X1>1))
//        L = L .* ( ((diag(P(X1, X2))./diag(Q(X1, X2))) .* (diag(P(X1, X1))./diag(Q(X1, X1))).^(K-1)) .* F + ~F )
//        F = X2<N+1
//        X1 = X2
//        i = i+1
//    end
//end
//
//disp(sum(Tn.*L)/nbSimulations)


























//
//// Definition de P (matrice de transition du systeme discretise)
//P = matriceTransition(lambda, mu)
//
//// Definition de Q (matrice de transition du changement de probabilite)
//Q = matriceTransition(new_lambda, new_mu)
//
//// Algorithme d'échantillonage préférentiel
n = 500
T = zeros(nbSimulations, 1)
F = ones(nbSimulations, 1)
X = ones(nbSimulations, 1)
while or(F)
    Y = grand(n, 'markov', Q, X)
    L = ones(nbSimulations, 1)
    for i=1:(n-1)
        A = diag(P(Y(:,i), Y(:,i+1))) ./ diag(Q(Y(:,i), Y(:, i+1)))
        L = L .* (A.* F + ~F)
//        for k=1:nbSimulations
//             L(k) = L(k) * (P(Y(k, i), Y(k, i+1))/Q(Y(k, i), Y(k, i+1)) * F(k) + ~F(k))
//        end
        T = T + h*F
        F = Y(:, i+1) < N+1
    end
    X = Y(:, $)
end

disp(sum(T.*L)/nbSimulations)











//
//esp=zeros(N,1);
//for m=1:m_max
//    for i=1:N
//        S=0;        
//        L = ones(M,1)
//        X = i*ones(M,1)
//        F = ones(M,1)
//        while or(F)
//            Xs = grand(1, 'markov', Q, X)
//            L = L.*(diag(P(X,Xs))./diag(Q(X,Xs)))
//            S = S + sum(F.*L)
//            X = Xs
//            F = X<(N+1)
//        end
//        esp(i) = max(1, S/M)
//    end
//    Q(1:N, :) = P(1:N, :)*diag(1+[esp;0])
//    Q = (diag(sum(Q, 'c'))^(-1))*Q
//    disp(m/m_max*100)
//end
//
//// Affichage
//disp(esp(1))
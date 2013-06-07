//------------------------------------------------------------------------------
// Changement de probabilité pour calculer E[Tn]
// 
//------------------------------------------------------------------------------

clear()
stacksize(150000000)

N = 40
nbSimulations = 100
h = 1

lambda=0.35
mu=0.4

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

function Y=f1(K)
    Y = (lambda*h*(1-lambda*h)^(K-1)) ./ (new_lambda*h*(1-new_lambda*h)^(K-1))
endfunction

function Y=f2(K)
    Y = (lambda*h*(1-(lambda+mu)*h)^(K-1)) ./ (new_lambda*h*(1-(new_lambda+new_mu)*h)^(K-1))
endfunction

function Y=f3(K)
    Y = (mu*h*(1-(lambda+mu)*h)^(K-1)) ./ (new_mu*h*(1-(new_lambda+new_mu)*h)^(K-1))
endfunction

stacksize(150000000)
n = 200
X = ones(nbSimulations, 1)
Tn = zeros(nbSimulations, 1)
F = ones(nbSimulations, 1)
L = ones(nbSimulations, 1)
Ones = ones(nbSimulations, 1)
m = zeros(nbSimulations, 1)
X2 = X
p = 1:nbSimulations
while or(F)
    i = 1
    T1 = grand(nbSimulations, n, 'bin', 1, new_lambda*h)
    T2 = grand(nbSimulations, n, 'bin', 1, (new_lambda+new_mu)*h)
    U = grand(nbSimulations, n, 'def')
    e = 1*(U<=new_lambda/(new_lambda+new_mu)) + (-1)*(U>new_lambda/(new_lambda+new_mu))
    while or(F) & (i<=n)
        K = (T1(:,i).*(X(:,$)==1) + T2(:,i).*(X(:,$)>1))
        Tn = Tn + h*F
//        X2 = X + F.*K(:, $).*((Ones.*(X==1)) + e(:,i).*(X>1))
        X = [X, X(:,$) + F.*K(:, $).*((Ones.*(X(:,$)==1)) + e(:,i).*(X(:,$)>1))]
//        for i=p(F)
//            L(i) = L(i) *Q(X(i), X2(i)) /  P(X(i), X2(i)) 
//        end
        //L = L .* ((f1(K).*(X==1) + (X>1).*(f2(K).*(e(:, i)==1) + f3(K).*(e(:, i)==-1))) .* F + ~F)
        F = X(:, $)<N+1
//        F = X2 < N+1
//        X = X2
        i = i+1
    end
    //disp('n')
end

L = ones(nbSimulations, 1)
for i=1:nbSimulations
    j = 1
    while X(i,j)<N+1
        L(i) = L(i) *Q(X(i, j), X(i, j+1)) /  P(X(i, j), X(i, j+1))
        j = j+1
    end
end

disp(sum(Tn.*L)/nbSimulations)

























// ----- /!\ INUTILE /!\ -----

//n = 100
//TnL = 0
//for m=1:nbSimulations
//   tn = 0
//   x1 = 1
//   x2 = 1
//   l = 1
//   while(x1 < N+1)
//       i = 1
//       T1 = grand(1, n, 'geom', new_lambda*h)
//       T2 = grand(1, n, 'geom', (new_lambda+new_mu)*h)
//       U = grand(1, n, 'def')
//       E = 1*(U<=new_lambda/(new_lambda+new_mu)) + (-1)*(U>new_lambda/(new_lambda+new_mu))
//       while(x1<N+1 & i<=n)
//           k = T1(i)*(x1==1) + T2(i)*(x1>1)
//           tn = tn + h*k
//           x2 = x1 + ((x1==1) + E(i)*(x1>1))
//           l = l * P(x1, x2)/Q(x1, x2) * (P(x1, x1)/Q(x2, x2))^(k-1)
//           i = i + 1
//       end
//   end 
//   TnL = [TnL, tn*l]
//end
//
//disp(sum(TnL)/nbSimulations)



























//
//// Definition de P (matrice de transition du systeme discretise)
//P = matriceTransition(lambda, mu)
//
//// Definition de Q (matrice de transition du changement de probabilite)
//Q = matriceTransition(new_lambda, new_mu)
//
//// Algorithme d'échantillonage préférentiel
//n = 500
//T = zeros(nbSimulations, 1)
//F = ones(nbSimulations, 1)
//X = ones(nbSimulations, 1)
//while or(F)
//    Y = grand(n, 'markov', Q, X)
//    L = ones(nbSimulations, 1)
//    for i=1:(n-1)
//        A = diag(P(Y(:,i), Y(:,i+1))) ./ diag(Q(Y(:,i), Y(:, i+1)))
//        L = L .* (A.* F + ~F)
////        for k=1:nbSimulations
////             L(k) = L(k) * (P(Y(k, i), Y(k, i+1))/Q(Y(k, i), Y(k, i+1)) * F(k) + ~F(k))
////        end
//        T = T + h*F
//        F = Y(:, i+1) < N+1
//    end
//    X = Y(:, $)
//end
//
//disp(sum(T.*L)/nbSimulations)











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
//------------------------------------------------------------------------------
// Changement de probabilité naïf en modifiant lambda et mu pour calculer E[Tn]
// 
//------------------------------------------------------------------------------

clear()
stacksize(150000000)

N = 3
nbSimulations = 1000
h = 1
n = 200

lambda=0.36
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

//function Y=f1(K)
//    Y = (lambda*h*(1-lambda*h)^(K-1)) ./ (new_lambda*h*(1-new_lambda*h)^(K-1))
//endfunction
//
//function Y=f2(K)
//    Y = (lambda*h*(1-(lambda+mu)*h)^(K-1)) ./ (new_lambda*h*(1-(new_lambda+new_mu)*h)^(K-1))
//endfunction
//
//function Y=f3(K)
//    Y = (mu*h*(1-(lambda+mu)*h)^(K-1)) ./ (new_mu*h*(1-(new_lambda+new_mu)*h)^(K-1))
//endfunction

stacksize(150000000)

//X = ones(nbSimulations, 1)
//Tn = zeros(nbSimulations, 1)
//F = ones(nbSimulations, 1)
//L = ones(nbSimulations, 1)
//Ones = ones(nbSimulations, 1)
//m = zeros(nbSimulations, 1)
//while or(F)
//    i = 1
//    T1 = grand(nbSimulations, n, 'bin', 1, new_lambda*h)
//    T2 = grand(nbSimulations, n, 'bin', 1, (new_lambda+new_mu)*h)
//    U = grand(nbSimulations, n, 'def')
//    e = 1*(U<=new_lambda/(new_lambda+new_mu)) + (-1)*(U>new_lambda/(new_lambda+new_mu))
//    while or(F) & (i<=n)
//        K = (T1(:,i).*(X(:,$)==1) + T2(:,i).*(X(:,$)>1))
//        Tn = Tn + h*F
//        X = [X, X(:,$) + F.*K(:, $).*((Ones.*(X(:,$)==1)) + e(:,i).*(X(:,$)>1))]
//        F = X(:, $)<N+1
//        i = i+1
//    end
//end
//
//for i=1:nbSimulations
//    j = 1
//    while X(i,j)<N+1
//        L(i) = L(i) *P(X(i, j), X(i, j+1)) /  Q(X(i, j), X(i, j+1))
//        j = j+1
//    end
//end

T = []
for j=1:nbSimulations
    X = 1
    Tn = 0
    L = 1
    while X<N+1
        i = 1
        T1 = h*grand(1, n, 'bin', 1, new_lambda*h)
        T2 = h*grand(1, n, 'bin', 1, (new_lambda+new_mu)*h)
        U =  grand(1, n, 'def')
        e = 1*(U<=new_lambda/(new_lambda+new_mu)) + (-1)*(U>new_lambda/(new_lambda+new_mu))
        while (X<N+1) & (i<=n)
            K = (T1(i).*(X==1) + T2(i).*(X>1))
            Tn = Tn + h
            X2 = X + K*((X==1) + e(i).*(X>1))
            L = L * P(X, X2) /  Q(X, X2)
            i = i+1
            X = X2
        end
    end
    T = [T, Tn*L]
end   

ValConf = 1.96*sqrt(variance(T, 'c'))/sqrt(nbSimulations)
E = sum(T, 'c')/nbSimulations
disp(E)
disp(E+ValConf, E-ValConf)
//------------------------------------------------------------------------------
// Implementation de l'algorithme du TD envoyé par F. Beynach-Georges dans notre
// cas
// 
// Rq: beaucoup plus précis mais pas de gain de temps
//------------------------------------------------------------------------------


clear;
stacksize(150000000);

h=0.1
lambda=0.45
mu=0.5
m_max=10;
N = 3;
M=10;

// Definition de P (matrice de transition)
P = zeros(N+1, N+1)
P(N+1, N+1)=1
P(1,1)=1-lambda*h
P(1,2)=lambda*h
for i=2:N
    P(i,i)=1-mu*h-lambda*h
    P(i,i-1)=mu*h
    P(i,i+1)=lambda*h
end

Q=P;

// Algorithme d'échantillonage préférentiel
esp=zeros(N,1);
for m=1:m_max
    for i=1:N
        S=0;        
        L = ones(M,1)
        X = i*ones(M,1)
        F = ones(M,1)
        while or(F)
            Xs = grand(1, 'markov', Q, X)
            L = L.*(diag(P(X,Xs))./diag(Q(X,Xs)))
            S = S + h*sum(F.*L)
            X = Xs
            F = X<(N+1)
        end
        esp(i) = max(1, S/M)
    end
    Q(1:N, :) = P(1:N, :)*diag(1+[esp;0])
    Q = (diag(sum(Q, 'c'))^(-1))*Q
    //disp(m/m_max*100)
end

// Affichage
disp(esp(1))










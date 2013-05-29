clear;
stacksize(150000000);

d=2;

cim=d+1;

eps=0.2;

alpha=2.2;

P=[0,eps,1-eps; 0, 1-eps^alpha, eps^alpha;  zeros(1,d),1];

mu=ones(d,1);


disp('Valeurs exactes :')

mu(1)=1-eps+eps*(1+eps^(-alpha));
mu(2)=eps^(-alpha);


for i=1:d
disp('  mu('+string(i)+') = '+string(mu(i)))
end

disp('Echantillonnage preferentiel :')

m_max=15;

M=25;

Q=P;



for m=1:m_max
	for i=1:d
	S=0;
		for l=1:M
			temps=0;
			X=i;
			L=1;
			while X<cim,
				Xs=grand(1,'markov', Q, X);
				L=L*P(X,Xs)/Q(X,Xs);
				temps=temps+L;
 				X=Xs;
 			end
 			S=S+temps;
 		end
 	mu(i)=max(1,S/M);
 	end
 	Q(1:d,:)=P(1:d,:)*diag(1+[mu;0]);
 	Q=diag((sum(Q,'c'))^(-1))*Q;
end



for i=1:d
disp('  mu('+string(i)+') = '+string(mu(i)))
end


disp('LGN classique :')


M2=M*m_max;


for i=1:d
	S=0;
		for l=1:M2
			temps=0;
			X=i;
			while X<cim,
				X=grand(1,'markov', P, X);
				temps=temps+1;
 			end
 			S=S+temps;
 		end
 	mu(i)=max(1,S/M2);
end


for i=1:d
disp('  mu('+string(i)+') = '+string(mu(i)))
end
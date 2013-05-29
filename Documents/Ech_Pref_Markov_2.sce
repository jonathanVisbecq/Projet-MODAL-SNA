clear;
stacksize(150000000);

d=19;

cim=d+1;

p=0.5;

D=[1,7,13,19];

l_D=length(D);

P=zeros(d+1,d+1);

for i=1:d
	for j=i:d
		P(i,j)=(1-p)^(j-i)*p;
	end
	P(i,$)=1-sum(P(i,1:d));
end

P(d+1,d+1)=1;

mu=ones(d,1);
mybeta=zeros(2,1);

disp('Valeurs exactes :')

mybeta(1)=(1+p*d)/(1-p);
mybeta(2)=-p/(1-p);


for i=1:2
disp('  beta('+string(i)+') = '+string(mybeta(i)))
end

disp('Echantillonnage preferentiel :')

E=[ones(d,1),(1:d)'];

E_D=E(D,:);

mutomybeta=inv(E_D'*E_D)*E_D';

m_max=2;

M=6;

Q=P;

for m=1:m_max
	for k=1:l_D
	S=0;
		for l=1:M
			temps=0;
			X=D(k);
			L=1;
			while X<cim,
				Xs=grand(1,'markov', Q, X);
				L=L*P(X,Xs)/Q(X,Xs);
				temps=temps+L;
 				X=Xs;
 			end
 			S=S+temps;
 		end
 	mu(D(k))=max(1,S/M);
 	end
 	mybeta=mutomybeta*(mu(D));
 	mu=max(1,(E*mybeta));
 	Q(1:d,:)=P(1:d,:)*diag(1+[mu;0]);
 	Q=diag((sum(Q,'c'))^(-1))*Q;
end



for i=1:2
disp('  beta('+string(i)+') = '+string(mybeta(i)))
end


disp('LGN classique :')


M2=M*m_max;


for k=1:l_D
	S=0;
		for l=1:M2
			temps=0;
			X=D(k);
			while X<cim,
				X=grand(1,'markov', P, X);
				temps=temps+1;
 			end
 			S=S+temps;
 		end
 	mu(D(k))=max(1,S/M2);
end
mybeta=mutomybeta*(mu(D));


for i=1:2
disp('  beta('+string(i)+') = '+string(mybeta(i)))
end
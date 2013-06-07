clear()
stacksize(150000000)

N = 5
nbSimulations = 1000
h = 1

lambda=0.55
mu=0.45

new_lambda = 0.6
new_mu = 0.4

n = 100

a = (lambda+mu)/(new_lambda+new_mu)
b = lambda + mu - new_lambda - new_mu

Tps = []
for m=1:nbSimulations
    X = 0
    L = 1
    Tn = 0
    while (X<N)
        j = 1
        T1 = grand(1, n, 'exp', 1/new_lambda)
        T2 = grand(1, n, 'exp', 1/(new_lambda+new_mu))
        U = grand(1, n, 'def')
        e = 1*(U<=new_lambda/(new_lambda+new_mu)) + (-1)*(U>new_lambda/(new_lambda+new_mu))
        while (X<N) & (j<n)
            T = (T1(j).*(X==0) + T2(j).*(X>0))
            Tn = Tn + T
            if X==0 then
                L = L * (lambda/new_lambda)*exp(-(lambda - new_lambda)*T)
            else
                L = L * ( (e(j)==1)*(lambda/new_lambda)/a + (e(j)==-1)*(mu/new_mu)/a )
            end        
            X = X + ((X==0) + (X>0)*e(j))
            j = j+1
        end
    end
    Tps = [Tps, Tn*L]
end

ValConf = 1.96*sqrt(variance(Tps))/sqrt(nbSimulations)
E = sum(Tps)/nbSimulations
disp(E)
disp(E+ValConf, E-ValConf)

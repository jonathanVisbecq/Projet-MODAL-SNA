//------------------------------------------------------------------------------
// Application de la transformee d'Esscher pour notre problème
//------------------------------------------------------------------------------

// Paramètres
lambda=0.45
mu=0.55
N = 10

// Nombre de simulations a effectuer
nbSimulations = 5000
// Nombre de variables aléatoires à simuler à chaque fois que nécessaire
n = 100

// Fonction f dans l'expression de la transformée d'Esscher
function y=f(x)
    if x==1 then
        y = 1/lambda
    elseif x==-1 then
        y = 1/mu
    else
        y = 0
    end
endfunction


function [p, e, eConf]=probaEspDepassementChgt(lambda, mu, N, nbSimulations, n)
    Ll = 0
    // Nouvelle intensité sous la transformation d'Esscher
    new_intensity = exp(f(1))*lambda + exp(f(-1))*mu
    // La nouvelle loi des saut est l'ancienne multipliée par un facteur exponentiel
    new_p = exp(f(1))*lambda/(exp(f(1))*lambda + exp(f(-1))*mu)
    
    Tps = []
    nb = 0
    for i=1:nbSimulations
        X = 1
        Xf = f(1)
        Tn = 0
        while (X>0) & (X<N)
            i = 1
            //T = h*grand(1, n, 'geom', (lambda+mu)*h)
            T = grand(1, n, 'exp', 1/(new_intensity))
            U =  grand(1, n, 'def')
            e = 1*(U<=new_p) + (-1)*(U>new_p)
            while (X>0) & (X<N) & (i<=n)
                Tn = Tn + T(i)
                X = X + e(i)
                Xf = Xf + f(e(i))
                i = i+1
            end
        end
        if X==N then
            nb = nb + 1
            //Tps = [Tps, Tn+h]
            // Calcul de la vraisemblance
            L = Xf - Tn*( lambda*(exp(f(1))-1) + mu*(exp(f(-1))-1) )
            L = exp(L)
            //disp(L)
            //disp(Tn)
            Tn = Tn/L
            Tps = [Tps, Tn]
        end
    end
    //disp(sum(Ll)/nbSimulations)
    p = nb/nbSimulations
    disp(length(Tps))
    e = sum(Tps)/length(Tps)
    eConf = 1.96*sqrt(variance(Tps))/sqrt(nbSimulations)
endfunction


[p, e, eConf]=probaEspDepassementChgt(lambda, mu, N, nbSimulations, n)
disp(e+eConf, e, e-eConf)
[p2, e2, eConf2] = probaEspDepassement(lambda, mu, N, nbSimulations, h, n)
disp(e2+eConf2, e2, e2-eConf2)



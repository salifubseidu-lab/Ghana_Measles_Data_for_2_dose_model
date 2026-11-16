function r = R0_fun(P)
%R0_FUN Basic reproduction number R0 = beta*eps/(k3*k40).
k = derived(P); r = P.beta*P.eps/(k.k3*k.k40);
end

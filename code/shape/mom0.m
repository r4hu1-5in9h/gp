function m0=mom0(mu,s2,lb,ub)

s=sqrt(s2);
arg_ub=(ub-mu)./s;
arg_lb=(lb-mu)./s;
m0=s.*(exp(logphi(arg_ub))-exp((logphi(arg_lb))));

end

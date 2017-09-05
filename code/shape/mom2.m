function m2=mom2(mu,s2,lb,ub)
f=normpdf([lb ub],mu,s2);
arg_ub=(ub-mu)./sqrt(2.*s2);
arg_lb=(lb-mu)./sqrt(2.*s2);
m2=0.5.*(s2+mu.^2).*(erf(arg_ub)-erf(arg_lb))...
    -s2.*((ub+mu).*f(2)-(lb+mu).*f(1));
end

function m1=mom1(mu,s2,lb,ub)
f=normpdf([lb ub],mu,s2);
arg_ub=(ub-mu)./sqrt(2.*s2);
arg_lb=(lb-mu)./sqrt(2.*s2);
m1=mu./2.*(erf(arg_ub)-erf(arg_lb))-s2.*(f(2)-f(1));
end


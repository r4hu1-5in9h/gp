function m0=mom0_old(mu,s2,lb,ub)
p=normcdf([lb ub],mu,s2);
m0=p(2)-p(1);
end

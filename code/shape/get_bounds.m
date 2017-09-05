function [lb ub] = get_bounds(w_price,w_inc)
%GET_BOUNDS obtains integration bounds for moment matching of
%   shape-constrained model

p=[1,(w_inc-1),w_price];
r=roots(p);

r = r(imag(r)==0); % Save only the real roots

lb=min(r); ub=max(r);

end


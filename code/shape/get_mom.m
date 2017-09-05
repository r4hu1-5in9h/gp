function [dot_Z, dot_mu, dot_s2]=get_mom(y,lb,ub,hyp,bar_mu,bar_s2)
%GET_MOM calculates the 0,1,2-moments in the shape-restricted model

s=exp(hyp.lik(1));
s2=s.^2;
tau=hyp.lik(2);

C=tau.*(1-tau)./s;
acute_C=exp(tau./s.*(bar_mu-y)+bar_s2./(2.*s2).*tau.^2);
grave_C=exp((tau-1)./s.*(bar_mu-y)+bar_s2./(2.*s2).*(tau-1).^2);

acute_mu=bar_mu+bar_s2./s.*tau;
grave_mu=acute_mu-bar_s2./s;

if y>=ub
    dot_Z=C.*acute_C.*mom0(acute_mu,bar_s2,lb,ub);
    dot_mu=C.*acute_C.*mom1(acute_mu,bar_s2,lb,ub)./dot_Z;
    dot_s2=C.*acute_C.*mom2(acute_mu,bar_s2,lb,ub)./dot_Z-dot_mu.^2;
elseif y<lb
    dot_Z=C.*grave_C.*mom0(grave_mu,bar_s2,lb,ub);
    dot_mu=C.*grave_C.*mom1(grave_mu,bar_s2,lb,ub)./dot_Z;
    dot_s2=C.*grave_C.*mom2(grave_mu,bar_s2,lb,ub)./dot_Z-dot_mu.^2;
else
    dot_Z=C.*acute_C.*mom0(acute_mu,bar_s2,lb,y)+...
        C.*grave_C.*mom0(grave_mu,bar_s2,y,ub);
    dot_mu=(C.*acute_C.*mom1(acute_mu,bar_s2,lb,y)+...
        C.*grave_C.*mom1(grave_mu,bar_s2,y,ub))...
        ./dot_Z;
    dot_s2=(C.*acute_C.*mom2(acute_mu,bar_s2,lb,y)+...
        C.*grave_C.*mom2(grave_mu,bar_s2,y,ub))...
        ./dot_Z-dot_mu.^2;
end

end
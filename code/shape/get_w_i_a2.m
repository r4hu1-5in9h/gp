function w_i_a=get_w_i_a2(X,alpha,hyp,i,a)

%GET_W_I_A returns the estimator for derivative of the GP with respect to
%   input component a evaluated at site i. The formula follows from
%   construction of the derivative process

[N,~]=size(X);
x_i=X(i,:)';

gamma_price2=exp(hyp.cov(1)).^2;
gamma_inc2=exp(hyp.cov(2)).^2;
gamma2=[gamma_price2, gamma_inc2]';
Gamma=diag(1./gamma2);
eta2=exp(hyp.cov(3)).^2;

% implementation2 : based on R&W prediction snippet in gp.m
if issparse(alpha)                  % handle things for sparse representations
    nz = alpha ~= 0;                                 % determine nonzero indices
else
    nz = true(size(alpha,1),1);
end               % non-sparse representation

K_a_xX=zeros(1,N);

for j=1:N
    x_j=X(j,:)';
    K_a_xX(j)=-eta2.*(x_i(a)-x_j(a))/gamma2(a).*exp(-0.5.*(x_i'*Gamma*x_j));
end

w_i_a=K_a_xX*full(alpha(nz,:));

%debug numerical problems
%w_i_a=min(w_i_a,1);
%w_i_a=max(w_i_a,0);

end


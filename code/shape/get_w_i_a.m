function w_i_a=get_w_i_a(X,K,ttau,tnu,hyp,i,a)

%GET_W_I_A returns the estimator for derivative of the GP with respect to
%   input component a evaluated at site i. The formula follows from
%   construction of the derivative process

[N,~]=size(X);
x_i=X(i,:)';
s=exp(hyp.lik(1));

gamma_price2=exp(hyp.cov(1)).^2;
gamma_inc2=exp(hyp.cov(2)).^2;
gamma2=[gamma_price2, gamma_inc2]';
Gamma=diag(1./gamma2);
eta2=exp(hyp.cov(3)).^2;

K_a_xX=zeros(1,N);

for j=1:N
    x_j=X(j,:)';
    K_a_xX(j)=-eta2.*(x_i(a)-x_j(a))/gamma2(a).*exp(-0.5.*(x_i'*Gamma*x_j));
end

% naive implementation
% w_i_a=(K_a_xX/(K_XX+tSigma))*tmu;

% see R&W algo 3.6, p59 of ch3
%sW = sqrt(ttau);                                     % compute Sigma and mu
%tS_half=diag(sW);
%L = chol(eye(N)+sW*sW'.*K);                         % L'*L=B=eye(n)+sW*K*sW
%z=tS_half*((L')\(L\(tS_half*K*tnu)));
%w_i_a=K_a_xX*(tnu-z);

% based on other code snippets in gpml
sW = sqrt(ttau);    

B=eye(N)+sW*sW'.*K;
%[BV,BD]=eig(B);
%disp(strcat('ttau: ',num2str(ttau)));
%disp(strcat('BD: ',num2str(BD)));
L = chol(B);                            % L'*L=B=eye(n)+sW*K*sW
%disp(strcat('L: ',num2str(L)));

alpha = tnu-sW.*solve_chol(L,sW.*(K*tnu));
w_i_a = K_a_xX*alpha;

% debug. this is where the issue is. also some NAs in tnu...

%debug numerical problems
%w_i_a=min(w_i_a,1);
%w_i_a=max(w_i_a,0);

end


clear;
startup;

%toy data
y=0.5; %note: this is y_i. suppress index throughout
i=3;
X=[1 2;7 5; 3 6; 9 0; 2 3];
mu=ones(5,1);
Sigma=eye(5);
hyp.lik=[.1 .5];
hyp.cov=[.05;.1;.08];
bar_mu=.5;
bar_s2=.5;
bar_tau=1./bar_s2;
bar_nu=bar_mu./bar_s2;

%step 1: obtain w_i_a forall i in {income, price}
w_price=get_w_i_a(X,mu,Sigma,hyp,i,1);
w_inc=get_w_i_a(X,mu,Sigma,hyp,i,2);


w_price=1;
w_inc=1;
%step 2: obtain q_i_lb, q_i_ub
[lb, ub]=get_bounds(w_price, w_inc);

%step 3: obtain dot_Z_i,dot_mu_i,dot_s2_i ==> dot_tau_i,dot_nu_i
[dot_Z, dot_mu, dot_s2]=get_mom(y,lb,ub,hyp,bar_mu,bar_s2);
dot_tau=1/dot_s2;
dot_nu=dot_mu/dot_s2;

%step 4: obtain ttau_i,teta_i
ttau=dot_tau-bar_tau;
tnu=dot_nu-bar_nu;
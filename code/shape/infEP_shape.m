function [post nlZ dnlZ] = infEP_shape(hyp, mean, cov, lik, x, y)

% Expectation Propagation approximation to the posterior Gaussian Process.
% The function takes a specified covariance function (see covFunctions.m) and
% likelihood function (see likFunctions.m), and is designed to be used with
% gp.m. See also infMethods.m. In the EP algorithm, the sites are
% updated in random order, for better performance when cases are ordered
% according to the targets.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2016-05-02.
%
% See also INFMETHODS.M.

persistent last_ttau last_tnu              % keep tilde parameters between calls
tol = 1e-4; max_sweep = 50; min_sweep = 2;     % tolerance to stop EP iterations

inf = 'infEP';
n = size(x,1);
if isnumeric(cov),  K = cov;                    % use provided covariance matrix
else [K,dK] = feval(cov{:},  hyp.cov,  x); end     % covariance matrix and deriv
if isnumeric(mean), m = mean;                         % use provided mean vector
else [m,dm] = feval(mean{:}, hyp.mean, x); end           % mean vector and deriv

% A note on naming: variables are given short but descriptive names in
% accordance with Rasmussen & Williams "GPs for Machine Learning" (2006): mu
% and s2 are mean and variance, nu and tau are natural parameters. A leading t
% means tilde, a subscript _ni means "not i" (for cavity parameters), or _n
% for a vector of cavity parameters. N(f|mu,Sigma) is the posterior. 
% Note the posterior is \tilde{r}(q)=N(q|\hat{mu},\hat{Sigma}) in
% Singh's notation

% marginal likelihood for ttau = tnu = zeros(n,1)
nlZ0 = -sum(feval(lik{:}, hyp.lik, y, m, diag(K), inf));
if any(size(last_ttau) ~= [n 1])      % find starting point for tilde parameters
    ttau = zeros(n,1); tnu  = zeros(n,1);        % init to zero if no better guess
    Sigma = K;                     % initialize Sigma and mu, the parameters of ..
    mu = m; nlZ = nlZ0;                  % .. the Gaussian posterior approximation
else
    ttau = last_ttau; tnu  = last_tnu;   % try the tilde values from previous call
    [Sigma,mu,L,alpha,nlZ] = epComputeParams(K, y, ttau, tnu, lik, hyp, m, inf);
    if nlZ > nlZ0                                           % if zero is better ..
        ttau = zeros(n,1); tnu  = zeros(n,1);       % .. then init with zero instead
        Sigma = K;                   % initialize Sigma and mu, the parameters of ..
        mu = m; nlZ = nlZ0;                % .. the Gaussian posterior approximation
    end
end

nlZ_old = Inf; sweep = 0;               % converged, max. sweeps or min. sweeps?
while (abs(nlZ-nlZ_old) > tol && sweep < max_sweep) || sweep<min_sweep
%while (sweep < max_sweep) || sweep<min_sweep    
    
    nlZ_old = nlZ; sweep = sweep+1;
    
    % added by Singh for moment matching
    %[Sigma,mu,L,alpha,nlZ] = epComputeParams(K, y, ttau, tnu, lik, hyp, m, inf);
    
    for i = randperm(n)       % iterate EP updates (in random order) over examples
        tau_ni = 1/Sigma(i,i)-ttau(i);      %  first find the cavity distribution ..
        nu_ni = mu(i)/Sigma(i,i)-tnu(i);                % .. params tau_ni and nu_ni
        
        % compute the desired derivatives of the indivdual log partition function
        
        % -----edit for shape restriction-----
        
        %step 0: keep old site param
        ttau_old = ttau(i); tnu_old = tnu(i);
        
        %step 1: obtain w_i_a forall i in {income, price}
        w_price=get_w_i_a(x,K,ttau,tnu,hyp,i,1);
        w_inc=get_w_i_a(x,K,ttau,tnu,hyp,i,2);
        %disp(w_price); %debug
        %disp(w_inc); %debug
        
        if or(isnan(w_price),isnan(w_inc))
            lb=[];
        else
            %step 2: obtain q_i_lb, q_i_ub
            [lb, ub]=get_bounds(w_price, w_inc);
        end

        
        if isempty(lb) %if no solution, skip this update
            ttau(i)=ttau_old;
            tnu(i)=tnu_old;
            disp('skip since no region satisfying constraint');
            
        else
            
            %disp(strcat('sweep: ',num2str(sweep)));
            %step 3: obtain dot_Z_i,dot_mu_i,dot_s2_i ==> dot_tau_i,dot_nu_i
            bar_mu=nu_ni/tau_ni;
            bar_s2=1/tau_ni;
            bar_tau=tau_ni;
            bar_nu=nu_ni;

            %disp(strcat('y= ',num2str(y(i))));
            %disp(strcat('lb= ',num2str(lb)));
            %disp(strcat('ub= ',num2str(ub)));
            %disp(strcat('bar_mu= ',num2str(bar_mu)));
            %disp(strcat('bar_s2= ',num2str(bar_s2)));
            %disp(strcat('hyp.lik(1)= ',num2str(hyp.lik(1))));
            %disp(strcat('hyp.lik(2)= ',num2str(hyp.lik(2))));
            
            [dot_Z, dot_mu, dot_s2]=get_mom(y(i),lb,ub,hyp,bar_mu,bar_s2);
           
            %disp(strcat('dot_Z= ',num2str(dot_Z)));    
            %disp(strcat('dot_mu= ',num2str(dot_mu)));
            %disp(strcat('dot_s2= ',num2str(dot_s2)));              
            
            if (dot_Z==0)
                
                ttau(i)=ttau_old;
                tnu(i)=tnu_old;
                disp('skip since dot_z=0');
                
            elseif (dot_s2==0)
                
                ttau(i)=ttau_old;
                tnu(i)=tnu_old;
                disp('skip since dot_s2=0');
            
            elseif (isnan(dot_s2))
                
                ttau(i)=0;
                tnu(i)=0;
                disp('reset since dot_s2=NaN');                
            
            else          
                dot_tau=1/dot_s2;
                dot_tau=max(dot_tau,0);
                dot_nu=dot_mu/dot_s2;

                %step 4: obtain ttau_i,tnu_i
                ttau(i)=dot_tau-bar_tau;
                ttau(i) = max(ttau(i),0);
                tnu(i)=dot_nu-bar_nu;
                
                %disp(strcat('bar_tau= ',num2str(bar_tau)));
                %disp(strcat('bar_nu= ',num2str(bar_nu)));

                disp(strcat('dot_s2= ',num2str(dot_s2)))
                %disp(strcat('dot_tau= ',num2str(dot_tau)));
                %disp(strcat('dot_nu= ',num2str(dot_nu)));
                
                %disp(strcat('ttau(i)= ',num2str(ttau(i))));
                %disp(strcat('tnu(i)= ',num2str(tnu(i))));
            
            end
            
            %smooth
            ttau(i)=(ttau(i)+ttau_old)/2;
            tnu(i)=(tnu(i)+tnu_old)/2;
            
        end
        
        % ------------------------------------

        
        dtt = ttau(i)-ttau_old; dtn = tnu(i)-tnu_old;      % rank-1 update Sigma ..
        si = Sigma(:,i); ci = dtt/(1+dtt*si(i));
        Sigma = Sigma - ci*si*si';                         % takes 70% of total time
        mu = mu - (ci*(mu(i)+si(i)*dtn)-dtn)*si;               % .. and recompute mu
    end
    % recompute since repeated rank-one updates can destroy numerical precision
    [Sigma,mu,L,alpha,nlZ] = epComputeParams(K, y, ttau, tnu, lik, hyp, m, inf);
end

if sweep == max_sweep && abs(nlZ-nlZ_old) > tol
    error('maximum number of sweeps exceeded in function infEP')
end

last_ttau = ttau; last_tnu = tnu;                       % remember for next call
post.alpha = alpha; post.sW = sqrt(ttau); post.L = L;  % return posterior params

if nargout>2                                           % do we want derivatives?
    dnlZ = hyp;                                   % allocate space for derivatives
    tau_n = 1./diag(Sigma)-ttau;             % compute the log marginal likelihood
    nu_n  = mu./diag(Sigma)-tnu;                    % vectors of cavity parameters
    sW = sqrt(ttau);
    F = alpha*alpha'-repmat(sW,1,n).*solve_chol(L,diag(sW));   % covariance hypers
    dnlZ.cov = -dK(F)/2;
    for i = 1:numel(hyp.lik)                                   % likelihood hypers
        dlik = feval(lik{:}, hyp.lik, y, nu_n./tau_n, 1./tau_n, inf, i);
        dnlZ.lik(i) = -sum(dlik);
    end
    [junk,dlZ] = feval(lik{:}, hyp.lik, y, nu_n./tau_n, 1./tau_n, inf);% mean hyps
    dnlZ.mean = -dm(dlZ);
end

% function to compute the parameters of the Gaussian approximation, Sigma and
% mu, and the negative log marginal likelihood, nlZ, from the current site
% parameters, ttau and tnu. Also returns L (useful for predictions).
function [Sigma,mu,L,alpha,nlZ] = epComputeParams(K,y,ttau,tnu,lik,hyp,m,inf)
n = length(y);                                      % number of training cases
sW = sqrt(ttau);                                        % compute Sigma and mu

%disp(strcat('sW: ',num2str(sW)));

B=eye(n)+sW*sW'.*K;
%[BV,BD]=eig(B);

%disp(strcat('BD: ',num2str(BD)));

L = chol(B);                            % L'*L=B=eye(n)+sW*K*sW

%disp(strcat('L: ',num2str(L)));

V = L'\(repmat(sW,1,n).*K);
Sigma = K - V'*V;
alpha = tnu-sW.*solve_chol(L,sW.*(K*tnu+m));
mu = K*alpha+m; v = diag(Sigma);

tau_n = 1./diag(Sigma)-ttau;             % compute the log marginal likelihood
nu_n  = mu./diag(Sigma)-tnu;                    % vectors of cavity parameters
lZ = feval(lik{:}, hyp.lik, y, nu_n./tau_n, 1./tau_n, inf);
p = tnu-m.*ttau; q = nu_n-m.*tau_n;                        % auxiliary vectors
nlZ = sum(log(diag(L))) - sum(lZ) - p'*Sigma*p/2 + (v'*p.^2)/2 ...
    - q'*((ttau./tau_n.*q-2*p).*v)/2 - sum(log(1+ttau./tau_n))/2;

function [varargout] = likALD(hyp, y, mu, s2, inf, i)

% likALD - Asymmetric Laplace density (ALD) likelihood.
% 
% likALD(hyp, y, mu, s2, inf, i)
% 
% The expression for the ALD is 
%   likALD(t) = q*(q-1)/sn * exp(-rho(y-t,q)/sn)
% where 
%   0<q<1 is the quantile of interest
%   y is the observation
%   sn is the scale
%   rho(x,q) = q*x         if x>=0
%              (q-1)*x     if x<0
%
% The hyperparameters are:
%   hyp = [  log(sn), q  ]
%
% In the Expectation-Propagation (EP) algorithm, this typically requires the mean
% and variance of the cavity field (context), here mu and s2 respectively.
% 
% inf specifies the type of inference (only 'infEP' supported).
% 
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likelihoods.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme. 
%
% Copyright (c) by Remi Barillec and Alexis Boukouvalas, 2012-04-25
% This is based on likLaplace, Copyright (c) by Carl Edward Rasmussen and Hannes
% Nickisch, 2010-07-21.
%
% See also likFunctions.m.

% LaTeX comments in the code indicates the corresponding terms in the technical
% report.

if nargin<2, varargout = {'2'}; return; end   % report number of hyperparameters

sn = exp(hyp(1));
q = hyp(2);

if nargin<5                              % prediction mode if inf is not present
  if numel(y)==0,  y = zeros(size(mu)); end
  s2zero = 1; if nargin>3, if norm(s2)>0, s2zero = 0; end, end         % s2==0 ?
  if s2zero                                         % log probability evaluation
    lp = -abs(y-mu)./b -log(2*b); s2 = 0;
  else                                                              % prediction
    lp = likALD(hyp, y, mu, s2, 'infEP');
  end
  ymu = {}; ys2 = {};
  if nargout>1
    ymu = mu;                                                   % first y moment
    if nargout>2
      ys2 = s2 + sn.^2;                                        % second y moment
    end
  end
  varargout = {lp,ymu,ys2};
  
else                                                            % inference mode
  switch inf 
  case 'infLaplace'
    error('Not implemented');
    
  case 'infEP'
    n = max([length(y),length(mu),length(s2),length(sn)]); on = ones(n,1);
    y = y.*on; mu = mu.*on; s2 = s2.*on; sn = sn.*on;             % vectors only
    fac = 1e3;          % factor between the widths of the two distributions ...
       % ... from when one considered a delta peak, we use 3 orders of magnitude
    idlik = fac*sn<sqrt(s2);                        % Likelihood is a delta peak - emission distribution is narrow
    idgau = fac*sqrt(s2)<sn;                          % Gaussian is a delta peak - posterior predictive from EP is narrow
    id = ~idgau & ~idlik;                          % interesting case in between
    
    if nargin<6                                             % no derivative mode
      lZ = zeros(n,1); dlZ = lZ; d2lZ = lZ;                    % allocate memory
      if any(idlik)
      error('check');
        [lZ(idlik),dlZ(idlik),d2lZ(idlik)] = ...
                                likGauss(log(s2(idlik))/2, mu(idlik), y(idlik)); % note wrong order y and mu - see likek Gauss
      end
      if any(idgau)
      error('check');
        [lZ(idgau),dlZ(idgau),d2lZ(idgau)] = ...
                                likALD(log(sn(idgau)), mu(idgau), y(idgau)); % note wrong order y and mu
      end
      if any(id)

        % substitution to obtain unit variance, zero mean Laplacian
        tmu = (mu(id)-y(id))./sn(id);     % \tilde{\mu}
        tvar = s2(id)./sn(id).^2;         % \tilde{\sigma}^2

        % log Z
        tstd = sqrt(tvar);          % \tilde{sigma}
        zm = tmu./tstd + q*tstd;    % \frac{m_{-}}{\tile{\sigme}}
        zp = zm - tstd;             % \frac{m_{+}}{\tile{\sigme}}
        
        Lm = tmu*q + 0.5*tvar*q^2;  
        Lp = - tmu - tvar * (q-0.5);
        
        am = logphi(-zm);
        ap = logphi( zp) - tmu - tvar * (q-0.5);

        lZ(id) = logsum2exp([ap,am]) + Lm - log(sn) + log(q*(1-q));
        
        if nargout>1
          % First derivative of log Z w.r.t \mu
          
          % log( N(z)/Phi(z) )
          lqm = -0.5*zm.^2 - 0.5*log(2*pi) - logphi(-zm);
          lqp = -0.5*zp.^2 - 0.5*log(2*pi) - logphi( zp);

          dam = - exp(lqm-0.5*log(s2(id)));
          dap =   exp(lqp-0.5*log(s2(id))) - 1./sn(id);
                    
          dLm = q./sn(id);
          
          % ( exp(ap).*dap + exp(am).*dam )./( exp(ap) + exp(am) )
          dlZ(id) = expABz_expAx([ap,am],[1;1],[dap,dam],[1;1]) + dLm;
          
          if nargout>2
            % Second derivative of log Z w.r.t \mu
            dapr = dap+1./sn(id);
            d2am = -zm ./ sqrt(s2(id)) * dam - dam.^2;
            d2ap = -zp ./ sqrt(s2(id)) * dapr - dapr.^2;
            
            bm = -zm/sqrt(s2(id)) * dam;            
            bp = -zp/sqrt(s2(id)) * (dap+1./sn(id)) - 2*dap./sn(id) - 1./sn(id)^2;
            
            d2lZ(id) = expABz_expAx([ap,am],[1;1],[bp,bm],[1;1]) - (dlZ(id)-dLm).^2;
          end
        end

      end
      varargout = {lZ,dlZ,d2lZ};
      
    else                                                       % derivative mode
      dlZhyp = zeros(n,1);
      if any(idlik)
        dlZhyp(idlik) = 0;
      end
      if any(idgau)
        error('Not implemented');
        dlZhyp(idgau) = ...
               likALD(log(sn(idgau)), mu(idgau), y(idgau), 'infLaplace', 1);
      end
      
      % Derivative of log(Z) w.r.t. to likelihood hyperparameters
      if any(id)
        if (i==1)
        % substitution to obtain unit variance, zero mean Laplacian
        tmu = (mu(id)-y(id))./sn(id);        
        tvar = s2(id)./sn(id).^2;
        tstd = sqrt(tvar);          % \tilde{sigma}
        sni = 1./sn(id);
        
        zm = tmu./tstd + q*tstd;     % \frac{m_{-}}{\tile{\sigme}}
        zp = zm - tstd;              % \frac{m_{+}}{\tile{\sigme}}
        dzm = -sni .* (q*tstd);
        dzp = dzm + sni .* tstd;
        
        dLm = - ( tmu*q + tvar * q^2 ) .* sni;
        dLp = tmu .* sni + 2*tvar .* sni * (q-0.5);
        
        % log( N(z)/Phi(z) )
        lqm = -0.5*zm.^2 - 0.5*log(2*pi) - logphi(-zm);
        lqp = -0.5*zp.^2 - 0.5*log(2*pi) - logphi( zp);
          
        am = logphi(-zm);
        ap = logphi( zp) - tmu - tvar * (q-0.5);
        dam = exp(lqm+log(-dzm));
        dap = exp(lqp+log(dzp)) + dLp;
        
        dlZhyp(id) = - sni + dLm + expABz_expAx([ap,am],[1;1],[dap,dam],[1;1]);
      else
        % Don't optimise w.r.t to the quantile - set gradient to zero
        dlZhyp(id) = zeros(n,1);
      end
      
      end
      varargout = {dlZhyp.*sn(id)};        % deriv. wrt hypers
    end

  case 'infVB'
    error('Not implemented')
  end
end

% computes y = log( sum(exp(x),2) ) in a numerically safe way by subtracting 
%  the row maximum to avoid cancelation after taking the exp
%  the sum is done along the rows
function [y,x] = logsum2exp(logx)
  N = size(logx,2);
  max_logx = max(logx,[],2);
  % we have all values in the log domain, and want to calculate a sum
  x = exp(logx-max_logx*ones(1,N));
  y = log(sum(x,2)) + max_logx;

% numerically safe implementation of f(t) = log(1-erf(t)) = log(erfc(t))
function f = lerfc(t)
  f  = zeros(size(t));
  tmin = 20; tmax = 25;
  ok = t<tmin;                               % log(1-erf(t)) is safe to evaluate
  bd = t>tmax;                                            % evaluate tight bound
  interp = ~ok & ~bd;                % linearly interpolate between both of them
  f(~ok) = log(2/sqrt(pi)) -t(~ok).^2 -log(t(~ok)+sqrt( t(~ok).^2+4/pi ));
  lam = 1./(1+exp( 12*(1/2-(t(interp)-tmin)/(tmax-tmin)) ));   % interp. weights
  f(interp) = lam.*f(interp) + (1-lam).*log(erfc( t(interp) ));
  f(ok) = f(ok) + log(erfc( t(ok) ));                                % safe eval

%  computes y = ( (exp(A).*B)*z ) ./ ( exp(A)*x ) in a numerically safe way
%  The function is not general in the sense that it yields correct values for
%  all types of inputs. We assume that the values are close together.
function y = expABz_expAx(A,x,B,z)
  N = size(A,2);  maxA = max(A,[],2);      % number of columns, max over columns
  A = A-maxA*ones(1,N);                                 % subtract maximum value
  y = ( (exp(A).*B)*z ) ./ ( exp(A)*x ); 

% safe implementation of the log of phi(x) = \int_{-\infty}^x N(f|0,1) df
% logphi(z) = log(normcdf(z))
function lp = logphi(z)
  lp = zeros(size(z));                                         % allocate memory
  zmin = -6.2; zmax = -5.5;
  ok = z>zmax;                                % safe evaluation for large values
  bd = z<zmin;                                                 % use asymptotics
  ip = ~ok & ~bd;                             % interpolate between both of them
  lam = 1./(1+exp( 25*(1/2-(z(ip)-zmin)/(zmax-zmin)) ));       % interp. weights
  lp( ok) = log( (1+erf(z(ok)/sqrt(2)))/2 );
  % use lower and upper bound acoording to Abramowitz&Stegun 7.1.13 for z<0
  % lower -log(pi)/2 -z.^2/2 -log( sqrt(z.^2/2+2   ) -z/sqrt(2) )
  % upper -log(pi)/2 -z.^2/2 -log( sqrt(z.^2/2+4/pi) -z/sqrt(2) )
  % the lower bound captures the asymptotics
  lp(~ok) = -log(pi)/2 -z(~ok).^2/2 -log( sqrt(z(~ok).^2/2+2)-z(~ok)/sqrt(2) );
  lp( ip) = (1-lam).*lp(ip) + lam.*log( (1+erf(z(ip)/sqrt(2)))/2 );

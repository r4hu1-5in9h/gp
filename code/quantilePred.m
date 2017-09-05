% Copyright (c) 2012 Alexis Boukouvalas, Remi Barillec, Dan Cornford.
% Disclaimer:
% This code comes without any warranty.
% The authors should be referenced in any publications quoting results from
% this code.
% Redistributions of this code must retain the above copyright notice and
% disclaimer.
% New versions of the quantile code will become available at
% http://wiki.aston.ac.uk/AlexisBoukouvalas
%
% To run this code, first install the GPML code library and run startup
% to make sure the path is setup correctly.covfunc = {'covMaterniso',5}; hyp2Original.cov = [0; 0;]; 
% The GPML code can be freely downloaded from
% http://www.gaussianprocess.org/gpml/code/matlab/doc/
%

%% Data
close all; 
clear all;
startup;

% Silverman's motorcycle benchmark dataset
data = [  2.4    0.0 ;  2.6   -1.3 ;  3.2   -2.7 ;  3.6    0.0 ;
          4.0   -2.7 ;  6.2   -2.7 ;  6.6   -2.7 ;  6.8   -1.3 ;
          7.8   -2.7 ;  8.2   -2.7 ;  8.8   -1.3 ;  8.8   -2.7 ;
          9.6   -2.7 ; 10.0   -2.7 ; 10.2   -5.4 ; 10.6   -2.7 ;
         11.0   -5.4 ; 11.4    0.0 ; 13.2   -2.7 ; 13.6   -2.7 ;
         13.8    0.0 ; 14.6  -13.3 ; 14.6   -5.4 ; 14.6   -5.4 ;
         14.6   -9.3 ; 14.6  -16.0 ; 14.6  -22.8 ; 14.8   -2.7 ;
         15.4  -22.8 ; 15.4  -32.1 ; 15.4  -53.5 ; 15.4  -54.9 ;
         15.6  -40.2 ; 15.6  -21.5 ; 15.8  -21.5 ; 15.8  -50.8 ;
         16.0  -42.9 ; 16.0  -26.8 ; 16.2  -21.5 ; 16.2  -50.8 ;
         16.2  -61.7 ; 16.4   -5.4 ; 16.4  -80.4 ; 16.6  -59.0 ;
         16.8  -71.0 ; 16.8  -91.1 ; 16.8  -77.7 ; 17.6  -37.5 ;
         17.6  -85.6 ; 17.6 -123.1 ; 17.6 -101.9 ; 17.8  -99.1 ;
         17.8 -104.4 ; 18.6 -112.5 ; 18.6  -50.8 ; 19.2 -123.1 ;
         19.4  -85.6 ; 19.4  -72.3 ; 19.6 -127.2 ; 20.2 -123.1 ;
         20.4 -117.9 ; 21.2 -134.0 ; 21.4 -101.9 ; 21.8 -108.4 ;
         22.0 -123.1 ; 23.2 -123.1 ; 23.4 -128.5 ; 24.0 -112.5 ;
         24.2  -95.1 ; 24.2  -81.8 ; 24.6  -53.5 ; 25.0  -64.4 ;
         25.0  -57.6 ; 25.4  -72.3 ; 25.4  -44.3 ; 25.6  -26.8 ;
         26.0   -5.4 ; 26.2 -107.1 ; 26.2  -21.5 ; 26.4  -65.6 ;
         27.0  -16.0 ; 27.2  -45.6 ; 27.2  -24.2 ; 27.2    9.5 ;
         27.6    4.0 ; 28.2   12.0 ; 28.4  -21.5 ; 28.4   37.5 ;
         28.6   46.9 ; 29.4  -17.4 ; 30.2   36.2 ; 31.0   75.0 ;
         31.2    8.1 ; 32.0   54.9 ; 32.0   48.2 ; 32.8   46.9 ;
         33.4   16.0 ; 33.8   45.6 ; 34.4    1.3 ; 34.8   75.0 ;
         35.2  -16.0 ; 35.2  -54.9 ; 35.4   69.6 ; 35.6   34.8 ;
         35.6   32.1 ; 36.2  -37.5 ; 36.2   22.8 ; 38.0   46.9 ;
         38.0   10.7 ; 39.2    5.4 ; 39.4   -1.3 ; 40.0  -21.5 ;
         40.4  -13.3 ; 41.6   30.8 ; 41.6  -10.7 ; 42.4   29.4 ;
         42.8    0.0 ; 42.8  -10.7 ; 43.0   14.7 ; 44.0   -1.3 ;
         44.4    0.0 ; 45.0   10.7 ; 46.6   10.7 ; 47.8  -26.8 ;
         47.8  -14.7 ; 48.8  -13.3 ; 50.6    0.0 ; 52.0   10.7 ;
         53.2  -14.7 ; 55.0   -2.7 ; 55.0   10.7 ; 55.4   -2.7 ;
         57.6   10.7 ]; 

     
     
XTrain = data(:,1);
YTrain = data(:,2);     
load('../data/data_BHP2.mat');

%% 1-dim X 
XTrain=log_y(1:200);
YTrain=log_q(1:200);

quantileToPredict = 0.75;
likfunc = @likALD; 
hyp.lik = [log(1) quantileToPredict];
covfunc = {'covSEard'}; 
hyp.cov = [0; 0;]; % Use sq exponential ARD covariance function

hyp = minimize(hyp, @gp, -10, @infEP, [], covfunc, likfunc, XTrain, YTrain); % maximum likelihood

% Prediction
%XTest = linspace(2,60,100)';
XTest = linspace(min(XTrain),max(XTrain),100)';
[ymu ys2 fmu fs2] = gp(hyp, @infEP, [], covfunc, likfunc, XTrain, YTrain, XTest);
% Only the noise-free predictions fmu,fs2 make sense in the quantile
% regression context since no noise likelihood is defined (see ICML 2012
% paper for more details).

% Plot
figure; hold on; title(sprintf('Quantile %g regression', quantileToPredict));
plot(XTrain, YTrain, 'om');            
plot(XTest, fmu, '-b'); 
plot(XTest, fmu + 2*sqrt(fs2),'--b'); 
plot(XTest, fmu - 2*sqrt(fs2),'--b'); 
        
n_obs=200;
x1=log_p(1:n_obs);
x2=log_y(1:n_obs);
y=log_q(1:n_obs);

XTrain=[x1 x2];
YTrain=y;


%% 2-dim X
load('../data/data_BHP2.mat');
n_obs=3640;
x1=log_p(1:n_obs);
x2=log_y(1:n_obs);
y=log_q(1:n_obs);

len=-10;

quantileToPredict = 0.25;
[x1_test,x2_test,fmu_block,fs2_block,hyp]=save_quantile(quantileToPredict,len,x1,x2,y);

quantileToPredict = 0.50;
[x1_test,x2_test,fmu_block,fs2_block,hyp]=save_quantile(quantileToPredict,len,x1,x2,y);

quantileToPredict = 0.75;
[x1_test,x2_test,fmu_block,fs2_block,hyp]=save_quantile(quantileToPredict,len,x1,x2,y);

%note: existing files have len=-10

%% see results
clear;
load('../data/results_0.5.mat');

CI_top=fmu_block+2.*sqrt(fs2_block);
CI_bottom=fmu_block-2.*sqrt(fs2_block);

scatter3(x1,x2,y,'filled');
hold on;
mesh(x1_test,x2_test,fmu_block);
mesh(x1_test,x2_test,CI_top);
mesh(x1_test,x2_test,CI_bottom);
xlabel('Price');
ylabel('Income');
zlabel('Demand');
view(30,10);
hold off;

%% work here

load('../data/data_BHP2.mat');

data=[log_q log_p log_y];
n_obs=100;
sample=datasample(data,n_obs,'Replace',false);
y=sample(:,1);
x1=sample(:,2);
x2=sample(:,3);

len=100;
quantiles = [0.25 0.50 0.75];
optimizers = {'sd','pr','lbfgs','sr1'};

for idz = 1:length(quantiles)
    
    quantileToPredict=quantiles(idz);
    
    tic
    %ticBytes(gcp);
    
    for idx = 1:length(optimizers)
        optimizer=optimizers{idx};
        [~,~,~,~,~]=save_quantile(quantileToPredict,len,optimizer,x1,x2,y);
    end
    
    %tocBytes(gcp)
    toc
    
end







tic
ticBytes(gcp);
for idx = 1:numel(quantiles)
    quantileToPredict = quantiles(idx);
    optimizer='sd';
    [~,~,~,~,~]=save_quantile(quantileToPredict,len,optimizer,x1,x2,y);
    optimizer='pr';
    [~,~,~,~,~]=save_quantile(quantileToPredict,len,optimizer,x1,x2,y);
    optimizer='lbfgs';
    [~,~,~,~,~]=save_quantile(quantileToPredict,len,optimizer,x1,x2,y);
    optimizer='sr1';
    [~,~,~,~,~]=save_quantile(quantileToPredict,len,optimizer,x1,x2,y);
end
tocBytes(gcp)
toc








load('../data/results_0.33_100_sr1.mat');
param_path=output.trace.param_path;
scatter3(param_path(1,:),param_path(2,:),param_path(3,:),'filled');
line(param_path(1,:),param_path(2,:),param_path(3,:));
%i reports # line searches until conv 

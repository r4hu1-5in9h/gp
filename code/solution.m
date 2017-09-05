%% Visualize tilted loss
clear;
loss_vis;

%% Visualize mean gasoline demand
prelim;


%% Compare optimization methods on subset of data
clear;
startup;
load('../data/data_BHP2.mat');
data=[share log_p log_y];

n_obs=100;
rng(3);
sample=datasample(data,n_obs,'Replace',false);
y=sample(:,1);
x1=sample(:,2);
x2=sample(:,3);

len=100;
quantiles = [0.25 0.50 0.75];
optimizers = {'sd','bb','pr','lbfgs'};
inf_string='infEP';

for idq = 1:length(quantiles)
    
    quantileToPredict=quantiles(idq);
    
    tic;
    
    figure;
    for ido = 1:length(optimizers)
        
        optimizer=optimizers{ido};
        fprintf('quantile: %3.2f optimizer: %s\n',quantileToPredict,optimizer);
        
        [~,~,~,~,~,output]=save_quantile(quantileToPredict,len,optimizer,x1,x2,y,inf_string);
        
        param_path=output.trace.param_path;
        n_iter=output.iterations;

        subplot(2,2,ido);
        scatter3(param_path(1,:),param_path(2,:),param_path(4,:),'filled');
        line(param_path(1,:),param_path(2,:),param_path(4,:));
        title(strcat(optimizers{ido},': ',num2str(n_iter)));
        xlabel('log(\gamma_{price})');
        ylabel('log(\gamma_{inc})');
        zlabel('log(\sigma)');
        view(30,30);
        
        [~,~,~,~,~,output]=save_quantile(quantileToPredict,len,optimizer,x1,x2,y,inf_string);
        
    end
    
    toc;
    
end


%% Implement chosen method on full data set

clear;
load('../data/data_BHP2.mat');
y=share;
x1=log_p;
x2=log_y;

len=10;
quantiles = [0.50]; %run one by one
optimizer='lbfgs';
inf_string='infEP';

%note: this takes many hours to run
train=false;

if train

    tic;
    %ticBytes(gcp);
    
    for ido = 1:numel(quantiles)
        quantileToPredict = quantiles(ido);
        [~,~,~,~,~,~]=save_quantile(quantileToPredict,len,optimizer,x1,x2,y,inf_string);
    end
    
    %tocBytes(gcp);
    toc;

end

%% Visualize quantiles of gasoline demand

clear;
load('../data/data_BHP2.mat');
y=log_q;
x1=log_p;
x2=log_y;

%debugging plots
load('../data/results_0.75_100_pr.mat');
plot_quantile3(x1,x2,y,x1_test,x2_test,fmu_block,fs2_block);
title('tau=0.75');

% convergence plot
fval=output.trace.fval;
plot(fval);

% log bayes factor
% fX is nlZ
% log(bayes factor of m over m')=logP(D|m)-logP(D|m')

load('../data/results_0.25_10_lbfgs.mat');
plot_quantile3(x1,x2,y,x1_test,x2_test,fmu_block,fs2_block);
%title('tau=0.25');

load('../data/results_0.5_10_lbfgs.mat');
plot_quantile3(x1,x2,y,x1_test,x2_test,fmu_block,fs2_block);
title('tau=0.50');

load('../data/results_0.75_10_lbfgs.mat');
plot_quantile3(x1,x2,y,x1_test,x2_test,fmu_block,fs2_block);
plot_quantile4(x1,x2,y,x1_test,x2_test,fmu_block,fs2_block);
plot_quantile5(x1,x2,y,x1_test,x2_test,fmu_block,fs2_block,47);

%%
% =============================
%  Subfunctions
% =============================
%%

%% loss_vis.m
% visualize different loss functions
%
% <include>loss_vis.m</include>
%
%%

%% prelim.m
% visualize mean gasoline demand
%
% <include>prelim.m</include>
%
%%

%% save_quantile.m
% train site parameters, tune hyperparameters, and save output
%
% <include>save_quantile.m</include>
%
%%

%% plot_quantile.m
% visualize quantile gasoline demand with credible interval
%
% <include>plot_quantile.m</include>
%
%%
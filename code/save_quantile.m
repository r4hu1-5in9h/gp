function [x1_test,x2_test,fmu_block,fs2_block,hyp,output]=save_quantile(quantileToPredict,len,optimizer,x1,x2,y,inf_string)
%SAVE_QUANTILE runs the main algorithm to train site parameters and tune
%hyperparameters. 

%It takes the following arguments:
%'quantileToPredict' quantile of interest
%'len' number of max permitted line searches
%'optimizer' string indicating which optimization algo to use
%'x1' vector of first components of training inputs
%'x2' vector of second components of training inputs
%'y' vector of training outputs
%'shape' indicates whether to impose shape restriction

%It saves the following objects in a .mat file in the data
%folder:
%'x1_test' vector of first components of the test inputs
%'x2_test' vector of second components of the test inputs
%'fmu_block' matrix of predicted latent quantile demand at test points
%'fs2_block' matrix of variances of those predictions
%'hyp' struct of tuned hyperparameters
%'fX' function value at solution
%'i' number of line searches
%'output'  struct of optimization information, including path of
%hyperparamter optimization

XTrain=[x1 x2];
YTrain=y;

if strcmp(inf_string,'infEP_shape')
    inf=@infEP_shape;
    shape_string='_restricted';
elseif strcmp(inf_string,'infEP_shape2')
    inf=@infEP_shape2;
    shape_string='_restricted2';
else
    inf=@infEP;
    shape_string='';
end

likfunc = @likALD; 
hyp.lik = [log(1) quantileToPredict];
covfunc = {'covSEard'}; 
hyp.cov = [0; 0; 0;]; % Use sq exponential ARD covariance function

if strcmp(optimizer,'wr')
    [hyp fX param_path i] = minimize_vis(hyp, @gp, len, inf, [], covfunc, likfunc, XTrain, YTrain); % maximum likelihood
    output = param_path;
elseif strcmp(optimizer,'sd')
    options=struct('MaxIter', len,'Method','sd');
    [hyp, fX, i, exitflag, output] = minimize_minfunc_vis(hyp, @gp, options, inf, [], covfunc, likfunc, XTrain, YTrain);
elseif strcmp(optimizer,'bb')
    options=struct('MaxIter', len,'Method','bb');
    [hyp, fX, i, exitflag, output] = minimize_minfunc_vis(hyp, @gp, options, inf, [], covfunc, likfunc, XTrain, YTrain);
elseif strcmp(optimizer,'fr')
    options=struct('MaxIter', len,'Method','pcg','cgUpdate',0);
    [hyp, fX, i, exitflag, output] = minimize_minfunc_vis(hyp, @gp, options, inf, [], covfunc, likfunc, XTrain, YTrain);
elseif strcmp(optimizer,'pr')
    options=struct('MaxIter', len,'Method','pcg');
    [hyp, fX, i, exitflag, output] = minimize_minfunc_vis(hyp, @gp, options, inf, [], covfunc, likfunc, XTrain, YTrain);
elseif strcmp(optimizer,'bfgs')
    options=struct('MaxIter', len,'Method','qnewton','qnUpdate',0);
    [hyp, fX, i, exitflag, output] = minimize_minfunc_vis(hyp, @gp, options, inf, [], covfunc, likfunc, XTrain, YTrain);
elseif strcmp(optimizer,'lbfgs')
    options=struct('MaxIter', len,'Method','lbfgs');
    [hyp, fX, i, exitflag, output] = minimize_minfunc_vis(hyp, @gp, options, inf, [], covfunc, likfunc, XTrain, YTrain);
elseif strcmp(optimizer,'sr1')
    options=struct('MaxIter', len,'Method','qnewton','qnUpdate',1);
    [hyp, fX, i, exitflag, output] = minimize_minfunc_vis(hyp, @gp, options, inf, [], covfunc, likfunc, XTrain, YTrain);
else
    fprintf('invalid optimizer \n');
end

n_test=100;
x1_test=linspace(min(x1),max(x1),n_test)';
x2_test=linspace(min(x2),max(x2),n_test)';

x1_stretch=repmat(x1_test,1,n_test)';
x1_stretch=x1_stretch(:);

x2_stretch=repmat(x2_test,n_test,1);

XTest=[x1_stretch x2_stretch];

[ymu, ys2, fmu, fs2,~, post] = gp(hyp, inf, [], covfunc, likfunc, XTrain, YTrain, XTest);

fmu_block=vec2mat(fmu,n_test)';
fs2_block=vec2mat(fs2,n_test)';

save(strcat('../data/results_',num2str(quantileToPredict,2),'_',num2str(len),'_',optimizer,shape_string,'.mat'), ...
    'x1_test','x2_test','fmu_block','fs2_block','post','hyp','fX','i','output')

end


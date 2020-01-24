%[train_data, train_outcome, train_list, test_data, test_outcome, test_list] = TCGA_readstats;
%% load data
clear all; close all;
load('/path/to/data/classification/val_data.mat')
load('/path/to/data/classification/classification\test_data.mat')
load('/path/to/data/classification/classification\train_data.mat')
load('/path/to/data/classification/classification\NIH_data.mat')


load('/path/to/data/classification/classification\val_outcome.mat')
load('/path/to/data/classification/classification\test_outcome.mat')
load('/path/to/data/classification/classification\train_outcome.mat')
load('/path/to/data/classification/classification\NIH_outcome.mat')

train_data = cat(1,train_data,val_data);
train_outcome = cat(1,train_outcome, val_outcome);


train_data = train_data(:,[1:49,140:229,500:587]);
test_data = test_data(:,[1:49,140:229,500:587]);
NIH_data = NIH_data(:,[1:49,140:229,500:587]);


% set hyperparameters we care about
numsplt = 3;
numlrncyc = 4000;
lrnrate = 0.001;
numpick = 20;
numvar = 15;
cost.ClassNames = [0; 1];
cost.ClassificationCosts = [0 1.382; 1 0];




%% feature selection

rng(1) % we always set random seed so we are always splitting the same way
t = templateTree('MaxNumSplits', numsplt, ...
    'surrogate','on', ...
    'NumVariablesToSample', 'all');
classificationEnsemble = fitcensemble(...
    train_data, ...
    train_outcome, ...
    'Method', 'AdaBoostM1', ...
    'NumLearningCycles', numlrncyc, ...
    'Learners', t, ...
    'LearnRate', lrnrate, ...
    'ClassNames', [0; 1], ...
    'KFold',10, ...
    'Prior','uniform', ...
    'Cost', cost);

allpicks = [];
for i = 1:numel(classificationEnsemble.Trained)    
    [imp,ma] = predictorImportance(classificationEnsemble.Trained{i});    
    %figure, bar(imp)
    [vals, picks] = maxk(imp,numpick);
    allpicks = cat(1,allpicks,picks);
end

[C,ia,ic] = unique(allpicks);
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];
[B,I] = sort(value_counts(:,2),'descend');
top_picks = find(B>1);

predictors = train_data(:,C(I(top_picks)));
%cd /dat
%% hyperparameter optimization - uncomment to test this
% params = hyperparameters('fitcensemble',train_data,train_outcome,'Tree');
% params(1,1).Optimize = false;    %do not optimize method
% params(2,1).Optimize = true;     %optimize num cycles
% params(2,1).Range = [200,1000];  %num cycles search range
% params(3,1).Optimize = true;     %optimize learning rate
% params(3,1).Range = [0.001,0.02]; %learning rate range
% params(4,1).Optimize = false;    %dont optimize mean leaf size for now
% params(5,1).Optimize = true;     %max num splits
% params(5,1).Range = [1,10];	 %range of num splits
% params(6,1).Optimize = false;	 %split criterion - dont change this
% params(7,1).Optimize = true;	 %num variables selected
% params(7,1).Range = [1,25];      %cap at 25
% 
% rng(1)
% t = templateTree('PredictorSelection','interaction-curvature');
% classificationEnsemble = fitcensemble(...
%     train_data, ...
%     train_outcome, ...
%     'Method', 'AdaBoostM1', ...
%     'OptimizeHyperparameters',params, ...
%     'ClassNames', [0; 1], ...
%     'Learners', t, ...
%     'Cost', cost);


%% 10 fold cross-val. we determine best model by OOB AUC
rng(1) % we always set random seed so we are always splitting the same way
t = templateTree('MaxNumSplits', numsplt, ...
    'PredictorSelection','interaction-curvature', ...
    'NumVariablesToSample',numvar);
classificationEnsemble = fitcensemble(...
    predictors, ...
    train_outcome, ...
    'Method', 'AdaBoostM1', ...
    'NumLearningCycles',numlrncyc, ...
    'Learners', t, ...
    'LearnRate', lrnrate, ...
    'ClassNames', [0; 1], ...
    'KFold',10, ...
    'Prior','uniform', ...
    'Cost', cost);

figure
plot(kfoldLoss(classificationEnsemble,'Mode','cumulative','LossFun','exponential'))
xlabel('Number of trees')
ylabel('Cross-validated exponential loss')

validationAccuracy = 1 - kfoldLoss(classificationEnsemble, 'LossFun', 'ClassifError')

[yFit,sFit] = kfoldPredict(classificationEnsemble);
figure, confusionchart(train_outcome,yFit);
[X,Y,T,trainAUC] = perfcurve(train_outcome,sFit(:,2),1);
trainAUC

%% once we are happy with our model, we apply it to entire training set for final feature selection
%  and then we apply that final model to the testing set to see how we did 

rng(1) % we always set random seed so we are always splitting the same way
t = templateTree('MaxNumSplits', numsplt, ...
    'PredictorSelection','interaction-curvature', ...
    'NumVariablesToSample',numvar);
finalModel = fitcensemble(...
    predictors, ...
    train_outcome, ...
    'Method', 'AdaBoostM1', ...
    'NumLearningCycles',numlrncyc, ...
    'Learners', t, ...
    'LearnRate', lrnrate, ...
    'ClassNames', [0; 1], ...
    'Prior','uniform', ...
    'Cost', cost);


[yfit,yscore] = predict(finalModel, test_data(:,C(I(top_picks))));
[X,Y,T,testAUC] = perfcurve(test_outcome,yscore(:,2),1);
testAUC
figure, y=confusionchart(test_outcome,yfit);
testAccuracy = (y.NormalizedValues(1)+y.NormalizedValues(4))/sum(y.NormalizedValues(:))
%%
[nfit,nscore] = predict(finalModel, NIH_data(:,C(I(top_picks))));
[X,Y,T,nihAUC] = perfcurve(NIH_outcome,nscore(:,2),1);
nihAUC
figure, y=confusionchart(NIH_outcome,nfit);
nihAccuracy = (y.NormalizedValues(1)+y.NormalizedValues(4))/sum(y.NormalizedValues(:))
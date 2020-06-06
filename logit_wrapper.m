%% Fit model to training set
[B,dev,stats] = mnrfit(Xtrain, Ytrain);

b = glmfit(Xtrain,Ytrain,'binomial','link','logit');

%% Predict probabilities for test set
Predicted_probabilities = glmval(b, Xtest,'logit');
%Predicted_probabilities = glmval(b,Xtest);

%% Count how many predictions were correct
for i=1:length(Ytest)
    predictedClassification(i) = categorical(round(Predicted_probabilities(i,1)+0.02));
end

%% real subject data
Predicted_probabilities_real = glmval(b, Xreal2,'logit');

for i=1:length(Xreal)
    predictedClassificationReal(i) = categorical(round(Predicted_probabilities_real(i,1)+0.02))
end


%% Confusion matrices

confusion = confusionmat(Ytest, predictedClassification)
TrueNeg = confusion(1,1);
FalsePos = confusion(1,2);
FalseNeg = confusion(2,1);
TruePos = confusion(2,2);
ActualNo = TrueNeg+FalsePos;
ActualYes = TruePos+FalseNeg;

accuracy = (TrueNeg + TruePos)/length(Ytest)
TruePosRate = TruePos/ActualYes;
TrueNegRate = TrueNeg/ActualNo;


    
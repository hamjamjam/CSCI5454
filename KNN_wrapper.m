k  = 11;
correct = 0;
clear predictedClassification;
M = ones(1,length(Xtrain(i,:)));

M = [1.0 0.8 2.0 0.8 0.7 9 2.2];

Xtrain_KNN = Xtrain;
Xtest_KNN = Xtest;
Xreal_KNN = Xreal;
for i = 1:length(M)
    M(i) = M(i)/mean(Xtrain(:,i));
    Xtrain_KNN(:,i) = Xtrain(:,i)*M(i);
    Xtest_KNN(:,i) = Xtest(:,i)*M(i);
    Xreal_KNN(:,i) = Xreal(:,i)*M(i);
end

for i = 1:length(Ytest)
    predictedClassification(i) = KNN(Xtrain_KNN, Ytrain, Xtest_KNN(i,:), k);
    if predictedClassification(i) == Ytest(i)
        correct = correct+1;
    end
end

for i = 1:length(Xreal)
    predictedClassificationReal(i) = KNN(Xtrain_KNN, Ytrain, Xreal_KNN(i,:),k);
end

confusion = confusionmat(Ytest, predictedClassification)
TrueNeg = confusion(1,1);
FalsePos = confusion(1,2);
FalseNeg = confusion(2,1);
TruePos = confusion(2,2);
ActualNo = FalseNeg+TruePos;
ActualYes = TrueNeg+FalsePos;

accuracy = (TrueNeg + TruePos)/length(Ytest)
TruePosRate = TruePos/ActualYes;
TrueNegRate = TrueNeg/ActualNo;

ind = ~isnan(Xtrain);
Xtrain2=Xtrain(ind);

    
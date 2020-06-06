
%% organize data
removal = 6;

Xall_dip = Xall;
Xall_dip(:,removal) = [];
Xreal_dip = Xreal;
Xreal_dip(:,removal) = [];


%% organize data
removal = 4;

Xall_dip(:,removal) = [];
Xreal_dip(:,removal) = [];

Mdl = fitcsvm(Xall_dip,Yall,'OutlierFraction',0.01);
cross = crossval(Mdl);

accuracy = 1-kfoldLoss(cross);

[label_real,score_real] = predict(Mdl,Xreal_dip);
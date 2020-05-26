%% Simulate a whole subject
%N = number of trials
%mu = underlying baseline mu
%sigma = underlying sigma
%levels = what SR noise levels are we operating on?
%withSR = do we have underlying SR? 0 for no SR 1 for SR
%initial_stim = initial stimulus of the task being simulated

function [thresholds, t, p]=simulateSubject(mu, sigma, withSR, N, initial_stim, ratios)

%% set up of threshold estimate simulations
guess_rate = 0.5; %at really small stim levels, 50% correct
lapse_rate = 0; %perfect/no lapses
plot_on = 0; %don't want any plots
sigma_guess = sigma; %starting point is real sigma

%We need to pick an underlying mu for each SR Noise level
if withSR == 1
    underlyingMus = mu*ratios; %this calls another function - for now, don't worry about that function
else
    underlyingMus = mu*ones(1,length(ratios));
end

%% initilaize list of thresholds
thresholds = ones(1,length(ratios));

%% simulate
%we are running that single threshold estimate simulation for every SR Noise Level
%in the 'levels' list because we want to assign a simulated threshold to to
%each SR Noise Level
for i = 1:length(ratios)
    mutemp = underlyingMus(i); %the mu we
    [X, Y, cor, lapses] = pest_mod_2int(N, mutemp, sigma, guess_rate, lapse_rate, plot_on, initial_stim);
    x = fminsearch(@(x) two_int_fit_simp(x, X, cor), [mutemp, sigma_guess]);
    if i == 1
        X_sham = X;
        cor_sham = cor;
        mu_sham = x(1);
    elseif i==2
        X_best = X;
        cor_best = cor;
        mu_best = x(1);
    else
        if x(1) < min(thresholds)
            X_best = X;
            cor_best = cor;
            mu_best = x(1);
        end
    end
    
    thresholds(i) = x(1); %add the simulated threshold to the list of thresholds
end

[t, p] = compareThresholds(X_sham,cor_sham,X_best,cor_best, mu_sham, mu_best, sigma_guess, N);

end
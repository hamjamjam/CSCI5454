clear all; close all; clc;
%
%% parameters that must match up to real life experiment
%these parameters are kept outside any for loops
%currently set for the visual task with VSR
N = 50;
levels = [0:0.1:1]; %VSR
initial_stim = 0.5; %visual task initially at 0.5 contrast

%% parameters that we change
%the rest of this file can go inside a for loop if we are simulating lots of
%subjects in one go
% withSR = zeros(1,2);
sim_number = 2000;
mumean = 0.15;
muSD = 0.05;

%% Setting up Features matrix
thresholds = zeros(length(levels), sim_number);
[m,n] = size(thresholds);
%  withSRi = zeros(1,sim_number);
%  withSRi(1:2:end) = 1;

underlyingMuRatios = underlyingMuRatiosVSR(levels);

for i = 1:sim_number
    running_sim_number = i
    features(i,1) = i;
    mu = normrnd(mumean, muSD); %pick mu
    if mu < 0
        mu = -mu;
    end
    withSR = randi([0 1],1,1); %withSRi(i); %we want to be able to choose if they have sr or not CHANGE TO 50/50
    sigma = mu*0.3; % we want this to be a function of mu for now 
    
    %% Simulate a subject and create features table
        %we get back a list of their thresholds for each stimulus level
        [thresholds(:,i),t,p]= simulateSubject(mu, sigma, withSR, N, initial_stim, underlyingMuRatios);
    %% Number of threshold values above sham, 2nd column of 'features'
        sham(i) = thresholds(1,i);
        features(i,2) = sum(thresholds(2:end,i) < sham(i));
    %% Number of direction changes ignoring first and last values, 3rd column of 'features'
        direction(:,i) = diff(thresholds(:,i));
        [m,n] = size(direction);
        directionchange = zeros(length(direction),2);
    
        for j = 1:(m-1)
            if (direction(j+1,i) > 0) && (direction(j,i) < 0) || (direction(j+1,i) < 0) && (direction(j,i) > 0)
                directionchange(j,1) = 1;
            else 
                directionchange(j,1) = 0;
            end
        end
        features(i,3) = sum(directionchange(:,1));
        
    %% Difference between lowest threshold value and sham
    
        %features(i,4) = sham(i) - min(thresholds(2:end,i));
        % this value will be positive unless sham is the lowest point, then
        % this value will be zero
    
    %% Ratio of best threshold to sham, 4th column of 'features'
        features(i,4) = min(thresholds(:,i))/sham(i);
        
    %% Avg of best and 2 neighboring, 5th column of 'features'
    
        if min(thresholds(:,i))==sham(i)
            minthresholds = [ thresholds(1,i), thresholds(2,i),  thresholds(3,i)];
        elseif min(thresholds(:,i))==thresholds(end,i)
            minthresholds = [ thresholds(end,i), thresholds(end-1,i),  thresholds(end-2,i)]; 
        else
            [m,n] = min(thresholds(1:end-1,i));
            minthresholds = [ thresholds(n,i), thresholds(n-1,i),  thresholds(n+1,i)];          
        end
        features(i,5) = mean(minthresholds)/sham(i); 
    %% Clarify if the subject has SR or no SR
       
        if withSR == 1
            SRvsNoSR(i,1) = 1;
        else
            SRvsNoSR(i,1) = 0;
        end

        %% Variance / STD of Threshold Values, 6th column of 'features'
        
        features(i,6) = sqrt(var(thresholds(:,i)))/sigma;
                
        features(i,7) = p; %p value of difference between best and sham
        
        features(i,8) = t; %p value of difference between best and sham

%     %% plot
%         %this is the SR plot for the subject we just simulated!
%         figure(1);
%         hold on
%         subplot(2, 5 ,(i)); % This needs to be changed depending on how many simulations are run
%         plot(levels, thresholds(:,i),'*--');
%         ylim([0 0.5]);
%         xticks(levels);
%         xlabel('SR Noise Level (mA)');
%         ylabel('Visual Contrast Threshold');
%         if withSR ==1
%             title('Simulated Subject with SR');
%         else
%             title('Simulated Subject without SR');
%         end
%         hold off
end
% %% Just the features as matrix and with vs without SR as vector
% Features = features(:,1:4);
% SRvNoSR = features(:,5);

%% Plotting with GPlotMatrix
% figure(2)
% A = features(:,[2 3 4 5]); % the columns 
% B = features(:,[6]); % the rows
% [~,ax] = gplotmatrix(A,B,SRvsNoSR(:)) % categorizing them based on if they have SR or not
% ax(1,1).YLabel.String = 'Standard Deviation of Thresholds';
% ax(1,1).XLabel.String = 'Num of Thresholds Above Sham';
% ax(1,2).XLabel.String = 'Num of Direction Changes';
% ax(1,3).XLabel.String = 'Ratio of Best to Sham';
% ax(1,4).XLabel.String = 'Avg of Best and Two Neighboring';

%% Create Labeled Table
%     %rowNames = {'A', 'B','c','d','e','f','g','h','i','j'};
%     colNames = {'Simulation Number', 'Number of Thresholds Above Sham', 'Number of Direction Changes', 'Diff Between Sham and Lowest','SR or No SR?'};
%     sTable = array2table(features,'VariableNames',colNames)

%% Split data into training and test sets
[row, col] = find(isnan(features));
features(row,:) = [];
SRvsNoSR(row) = [];
split_ind = floor(0.8   *length(features)); %TEMPORARILY CHANGING TO 50/50 FOR CONFUSION MATRIX

Xtrain = features(1:split_ind,2:end);
Ytrain = categorical(SRvsNoSR(1:split_ind,1));
Xtest = features((split_ind+1):end,2:end);
Ytest = categorical(SRvsNoSR((split_ind+1):end,1));

Xall = features(:,2:end);
Yall = categorical(SRvsNoSR(:,1));
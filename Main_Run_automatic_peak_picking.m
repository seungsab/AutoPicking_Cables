%% Fully-automated peak-picking method using modified-AMPD and MAD after baseline-correction
% The proposed method composed of the following three methods
% M-AMPD is the modified AMPD by padding to mitigate edge effect in original AMPD.
% MAD: median absolute deviation for robust outlier detection in presence of oultiers
% Baseline-correction: Baseline estimation algorithm in Raman Spectra
% CODED BY S.S.JIN (2020.08.28)

%% INITIALIZE
clear all; clc; close all;

%% Add Path 
addpath(genpath([pwd '\Utils']));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select the class
Class_Type = 'Easy';
% Class_Type = 'Intermediate';
% Class_Type = 'Hard';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load PSDs from the selected class
load([pwd '\Data\Faster_RCNN_Box_final.mat']);

switch Class_Type
    case 'Easy'
        fn_expert = 'Peak_expert_Easy.mat';
        fn_f_PSD = 'f_PSD_data_class_easy_data.mat';
        fn_Faster_RCNN_time = 'time_consumption_easy.txt';
        Final_peak_F_RCNN = easy_2;
        
    case 'Intermediate'
        fn_expert = 'Peak_expert_Intermediate.mat';
        fn_f_PSD = 'f_PSD_data_class_intermediate_data.mat';
        fn_Faster_RCNN_time = 'time_consumption_inter.txt';
        Final_peak_F_RCNN = intermediate_2;
    
    case 'Hard'
        fn_expert = 'Peak_expert_Hard.mat';
        fn_f_PSD = 'f_PSD_data_class_hard_data.mat';
        fn_Faster_RCNN_time = 'time_consumption_hard.txt';
        Final_peak_F_RCNN = hard_2;
end

% Set sampling frquency
fs = 100;

%% Import results from experts and Faster R-CNN
load([pwd '\Data\' fn_expert]); load([pwd '\Data\' fn_f_PSD]);

% Computing time for prediction via Faster R-CNN
computing_time_RCNN= load([pwd '\Data\' fn_Faster_RCNN_time]);

%% Perform Comparative Study for Peak-picking
Peak = {}; computing_time = [];
N_peak_expert = []; Accuracy = [];

for i_time = 1:size(PSD,1)
    legend_str = {}; n_ind = 1; Nind_COM = 1;
    f= Freq{i_time,1}; Cxx = PSD{i_time,1};
    
    %% Plot the PSD graph
    figure(1); cla
    plot(f,Cxx,'-','linewidth',2.5); % PSD
    xlabel('Frequency (Hz)','fontsize',15,'fontweight','bold');
    ylabel('PSD (in log scale)','fontsize',15,'fontweight','bold')
    title(['Class: ' Class_Type '  //  Sample Index: ' num2str(i_time)],'fontsize',15,'fontweight','bold')
    set(gca,'fontsize',25,'fontweight','bold')
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.8]);
    hold on;
    
    
    %% Method #1: Peak-picking via Faster R-CNN (results from Tensorflow)
    % Import the results from tensorflow
    Peak{i_time,n_ind} = Final_peak_F_RCNN{i_time,1};
    f_FRCNN = Peak{i_time,n_ind}(:,1);
    
    % Obtain Computing time
    computing_time(i_time,Nind_COM) = computing_time_RCNN(i_time); Nind_COM = Nind_COM + 1;    
    
    
    % Plot - Faster R-CNN
    ind_obj(n_ind) = plot(Peak{i_time,n_ind}(:,1),Peak{i_time,n_ind}(:,2),'gs',...
        'linewidth',3,'markersize',30); % Peaks selected
    legend_str = [legend_str {' Faster R-CNN'}];    
    n_ind = n_ind + 1;
        
    %% Method #2: Modifed Automatic Multiscale-based Peak Detection (M-AMPD)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run Modified Automatic Multiscale Peak Detection (AMPD)
    tStart = tic;
    ind_Mampd=Modified_ampd(Cxx,{fn_f_PSD, f}); 
    Cxx_ampd= Cxx(ind_Mampd); f_peaks_ampd = f(ind_Mampd);
    
    % Obtain Computing time
    computing_time(i_time,Nind_COM) = toc(tStart); Nind_COM = Nind_COM + 1;
    
    Peak{i_time,n_ind} = [f_peaks_ampd Cxx_ampd];
    f_AMPD = Peak{i_time,n_ind}(:,1);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Plot - Modifed AMPD
    ind_obj(n_ind) = plot(Peak{i_time,n_ind}(:,1),Peak{i_time,n_ind}(:,2),'ko',...
        'markersize',18,'linewidth',3); % Peaks selected
    legend_str = [legend_str {' Modified AMPD'}];
    n_ind = n_ind + 1;
    
    %% Method #3:  MAD after Baseline Removal
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tStart = tic;
    % Run Baseline-correction Algorithm 
    [Base, Cxx1]=baseline(Cxx,{fn_f_PSD, f});
    
    % Run Median Absolute Deviation (MAD) for Baseline removal (Cxx1)
    C = median(Cxx1); % Median
    c=-1/(sqrt(2)*erfcinv(3/2)); % Scaling factor for outlier detection using MAD
    threshold(1) = C+3*c*median(abs(Cxx1-median(Cxx1))); % Upper Threshold for outlier classification
    threshold(2) = C-3*c*median(abs(Cxx1-median(Cxx1))); % Lower Threshold for outlier classification
    
    ind_MAD1 = find(Cxx1>threshold(1)); ind_MAD2 = find(Cxx1<threshold(2));
    ind_MAD = union(ind_MAD1,ind_MAD2);
    
    temp_Cxx = min(Cxx)*ones(size(Cxx1,1),1);
    temp_Cxx(ind_MAD) = Cxx(ind_MAD);
    [~,ind_MAD]=findpeaks(temp_Cxx); % Find the maxima by find peaks
    
    Cxx_MAD_BS = Cxx(ind_MAD); f_peaks_MAD_BS = f(ind_MAD);
    
    % Obtain Computing time
    computing_time(i_time,Nind_COM) = toc(tStart); Nind_COM = Nind_COM + 1;
    
    Peak{i_time,n_ind} = [f_peaks_MAD_BS Cxx_MAD_BS];
    f_MAD_BS = Peak{i_time,n_ind}(:,1);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
    % Plot - MAD After BS
    ind_obj(n_ind) = plot(Peak{i_time,n_ind}(:,1),Peak{i_time,n_ind}(:,2),'rx',...
        'markersize',15,'linewidth',3); % Peaks selected
    legend_str = [legend_str {' MAD after Baseline-correction'}];
    n_ind = n_ind + 1;

    %% Method #4 Proposed method: Peak-picking via AMPD + MAD after Baseline Removal
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) Run Modified Automatic Multiscale Peak Detection (AMPD)
    ind_Mampd=Modified_ampd(Cxx,{fn_f_PSD, f});
    
    % 2) Run Baseline-correction Algorithm
    [Base, Cxx1]=baseline(Cxx,{fn_f_PSD, f});
    
    % 3) Run Median Absolute Deviation (MAD) for Baseline removal (Cxx1)
    C = median(Cxx1); % Median
    c=-1/(sqrt(2)*erfcinv(3/2)); % Scaling factor for outlier detection using MAD
    threshold(1) = C+3*c*median(abs(Cxx1-median(Cxx1))); % Upper Threshold for outlier classification
    threshold(2) = C-3*c*median(abs(Cxx1-median(Cxx1))); % Lower Threshold for outlier classification
    
    ind_MAD1 = find(Cxx1>threshold(1)); ind_MAD2 = find(Cxx1<threshold(2));
    ind_MAD = union(ind_MAD1,ind_MAD2);
    
    temp_Cxx = min(Cxx)*ones(size(Cxx,1),1);
    temp_Cxx(ind_MAD) = Cxx(ind_MAD);
    [~,ind_MAD]=findpeaks(temp_Cxx); % Find the maxima by find peaks
    
    % 4) Information integration 
    ind3 = intersect(ind_Mampd,ind_MAD); % peaks selected by both methods
    Cxx_peaks1 = Cxx(ind3); f_peaks1 = f(ind3);
        
    % Obtain Computing time by simply adding computing times in Method #2
    % (M-AMPD) and Method #3 (MAD after Baseline-correction)
    computing_time(i_time,Nind_COM) = sum(computing_time(i_time,2:3));
    Nind_COM = Nind_COM + 1;
    
    Peak{i_time,n_ind} = [f_peaks1 Cxx_peaks1];
    f_AMPD_MAD = Peak{i_time,n_ind}(:,1);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Plot - Proposed mnethod
    ind_obj(n_ind) = plot(Peak{i_time,n_ind}(:,1),Peak{i_time,n_ind}(:,2),'bo',...
        'markersize',40,'linewidth',3); % Peaks selected
    legend_str = [legend_str {' Proposed Method'}];
    n_ind = n_ind + 1;
    
    
    %% Plot Peaks selected by experts as Reference
    Peak{i_time,n_ind} = Final_peak_expert{i_time,1};
    f_exeprt = Peak{i_time,n_ind}(:,1);
    
    YLIM = get(gca,'ylim'); if YLIM(2)>0, YLIM(2)=YLIM(2)*1.3; else, YLIM(2)=YLIM(2)*0.7; end
    for jj = 1:size(Final_peak_expert{i_time,1},1)
        ind_obj(n_ind) = plot(Peak{i_time,n_ind}(jj,1)*[1 1],YLIM,':','linewidth',2,...
            'color',[0.85 0.33 0.1]);
    end
    legend_str = [legend_str {' Expert''s Perception'}];
    n_ind = n_ind + 1;
    
    %% SAVE FIGURE
    leg = legend(ind_obj,legend_str,'fontsize',20,'fontweight','bold');
    set(leg, 'Location', 'North'); set(leg, 'Location', 'Best'); axis tight
    xlim([0 12.5]);
    
    %% COMPUTE ACCURACY
    % # of True Peaks Selected by Experts
    N_peak_expert(i_time,1) = length(f_exeprt);
    
    % Due to different floating points in Python (Faster R-CNN) and MATLAB (Others)
    % Rounding the frequency to remove this influence
    f_FRCNN1 = round(f_FRCNN,4,'significant');
    f_AMPD1 = round(f_AMPD,4,'significant');
    f_MAD_BS1 = round(f_MAD_BS,4,'significant');
    f_AMPD_MAD_BS1 = round(f_AMPD_MAD,4,'significant');
    f_exeprt1 = round(f_exeprt,4,'significant');
    
    % Compute # True Peak Detected (TP)
    Accuracy.TP(i_time,1) = size(intersect(f_exeprt1,f_FRCNN1),1); % Between expert and Faster R-CNN
    Accuracy.TP(i_time,2) = size(intersect(f_exeprt1,f_AMPD1),1); % Between expert and M-AMPD
    Accuracy.TP(i_time,3) = size(intersect(f_exeprt1,f_MAD_BS1),1); % Between expert and MAD_BS
    Accuracy.TP(i_time,4) = size(intersect(f_exeprt1,f_AMPD_MAD_BS1),1); % Between expert and (AMPD + MAD_BS)
    
    % Compute # True Peak Missing (TN)
    Accuracy.TN(i_time,1) = size(setdiff(f_exeprt1,f_FRCNN1),1); % Between expert and Faster R-CNN
    Accuracy.TN(i_time,2) = size(setdiff(f_exeprt1,f_AMPD1),1); % Between expert and M-AMPD
    Accuracy.TN(i_time,3) = size(setdiff(f_exeprt1,f_MAD_BS1),1); % Between expert and MAD_BS
    Accuracy.TN(i_time,4) = size(setdiff(f_exeprt1,f_AMPD_MAD_BS1),1); % Between expert and (AMPD + MAD_BS)
    
    % Compute # False Peak Detected (FP)
    Accuracy.FP(i_time,1) = size(setdiff(f_FRCNN1,f_exeprt1),1); % Between expert and Faster R-CNN
    Accuracy.FP(i_time,2) = size(setdiff(f_AMPD1,f_exeprt1),1); % Between expert and M-AMPD
    Accuracy.FP(i_time,3) = size(setdiff(f_MAD_BS1,f_exeprt1),1); % Between expert and MAD_BS
    Accuracy.FP(i_time,4) = size(setdiff(f_AMPD_MAD_BS1,f_exeprt1),1); % Between expert and (AMPD + MAD_BS)    
end
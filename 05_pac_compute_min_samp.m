%% Calculate the minimum sample available per participant across conditions before computing the PAC;
clearvars;clc; 
clean_dir = 'XXX';
pac_out = 'XXX';
peak_dir = 'XXX';
pac_tmp = 'XXX';

load stats_beta.mat; load all_labels.mat;

subj_list_clean = dir([clean_dir,'*.mat']);
subj_list_randtrl = dir([pac_out,'rand_trials\','*.mat']);


% cycle over participants;
count_subj = 0; group_dist = [];

for s = 1 : length(subj_list_clean)/2

count_subj = count_subj+1;

% radnomize rand at each iteration/subj;
rng('shuffle'); 

% cycle over tasks;
expe = {'L1N','L1S'};
count = 0; idx_name = [];
for item = 1 : length(expe)
    task = expe{item};
    idx = [];
    for ii = 1:length(subj_list_clean)
        idx(1,ii) = strcmp(subj_list_clean(s+count).name,subj_list_clean(ii).name);
    end
     idx_name(1,item) = find(idx==1);
     count = count+18;   
end
    
    % load data_clean and peak info of the participant;
    name1 = subj_list_clean(idx_name(1,1)).name;    
    name3 = subj_list_randtrl(idx_name(1,1)).name;
    load([clean_dir,name1]); load([pac_out,'rand_trials\',name3]); 
    data_N = data_clean; rand_trials_N = unique([rand_trial_nc;rand_trial_c]);
    name1 = subj_list_clean(idx_name(1,2)).name;    
    name3 = subj_list_randtrl(idx_name(1,2)).name;
    load([clean_dir,name1]); load([pac_out,'rand_trials\',name3]);   
    data_S = data_clean; rand_trials_S = unique([rand_trial_nc;rand_trial_c]);

    % keep only hit trials;
    cfg = [];
    cfg.channel = 'all';
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.hit,1), data_N.trialinfo, 'UniformOutput', false)));
    data_N = ft_selectdata(cfg,data_N);
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.hit,1), data_S.trialinfo, 'UniformOutput', false)));
    data_S = ft_selectdata(cfg,data_S);    

    % now subsample the trials from Cognate and Non-cognate conditions;  
    min_samp = round(min([length(data_N.trialinfo) length(data_S.trialinfo)]));     
    idx_n = 1:length(data_N.trialinfo); 
    idx_s = 1:length(data_S.trialinfo); 
    rand_trial_n = sortrows(idx_n(randperm(length(idx_n),min_samp)));
    rand_trial_s = sortrows(idx_s(randperm(length(idx_s),min_samp)));
    %select subsamples trials;
    trial_n_temp = data_N.trial(1,rand_trial_n);
    trial_s_temp = data_S.trial(1,rand_trial_s);
    %corresponding trial inf;  
    info_n_temp = data_N.trialinfo(rand_trial_n);
    for i = 1 : length(info_n_temp)
        info_n_temp{i,1}.task(1,1) = 1;
    end
    info_s_temp = data_S.trialinfo(rand_trial_s); 
    for i = 1 : length(info_s_temp)
        info_s_temp{i,1}.task(1,1) = 2;
    end
    % create temporary structure;
    data_temp = [];
    data_temp.trial = cat(2,trial_n_temp,trial_s_temp);
    data_temp.trialinfo = cat(1,info_n_temp,info_s_temp);
    data_temp.fsample = data_N.fsample;
    data_temp.label = data_N.label;
    % add corresponding time info;
    time_temp = data_N.time{1,1};
    for tm = 1: length(data_temp.trialinfo)     
        data_temp.time{1,tm} = time_temp;
    end

    cfg = [];
    cfg.channel = [stats_beta.L1N.roi;stats_beta.L1S.roi];
    data_temp = ft_selectdata(cfg,data_temp);

    % get minimum number of trials available between cognates (cond = 1) and non cognate conditions (cond = 2);
    trlinfo = zeros(numel(data_temp.trialinfo),1);
    for j = 1 : numel(trlinfo)
        trlinfo(j,1) = j;
        trlinfo(j,2) = data_temp.trialinfo{j}.trigger;
        trlinfo(j,3) = data_temp.trialinfo{j}.condition;
        trlinfo(j,4) = data_temp.trialinfo{j}.task;
    end

    % balance the number of trials in each conditions (min trial available across the 4 conditions);
    min_samp = round(min([sum(trlinfo(:,3)==1 &trlinfo(:,4)==1) sum(trlinfo(:,3)==2 &trlinfo(:,4)==1) ...
                          sum(trlinfo(:,3)==1 &trlinfo(:,4)==2) sum(trlinfo(:,3)==2 &trlinfo(:,4)==2)])); 


    allsub_minsample(s,1) = min_samp;


end

save([pac_out,'allsub_minsample'],'allsub_minsample');

%%
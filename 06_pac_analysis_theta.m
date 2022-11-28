% PAC Analysis;
clearvars;clc;
clean_dir = 'XXX';
tfr_dir = 'XXX';
pac_out = 'XXX';
pac_out2 = 'XXX';
peak_dir = 'XXX';

load stats_beta.mat;

subj_list = dir([clean_dir,'*.mat']);
subj_list_peak = dir([peak_dir,'peak_infos\','*.mat']);

load([pac_out,'allsub_minsample.mat']);
allsub_minsample = repmat(allsub_minsample,2,1);

% loop over tasks;
expe = {'L1N','L1S'};
close all;

for item = 1 : length(expe)

% select file names only from expe(item);
task = expe{item};

clear idx
for ii = 1 : length(subj_list)
    idx(1,ii) = contains(subj_list(ii).name,task);
end
idx = find(idx==1);

% loop over participants;
for s = idx 
  
% radnomize rand at each iteration/subj;
rng('shuffle'); 

% load data_clean and peak info of the participant;
name1 = subj_list(s).name;    
name2 = subj_list_peak(s).name;
load([clean_dir,name1]); load([peak_dir,'peak_infos\',name2]);

% keep trials and channels of interest;
cfg = [];
cfg.channel = 'all';
cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.hit,1), data_clean.trialinfo, 'UniformOutput', false)));
data_hits = ft_selectdata(cfg,data_clean);
    
% get minimum number of trials available between cognates (cond = 1) and non cognate conditions (cond = 2);
trlinfo = zeros(numel(data_hits.trialinfo),1);
for j = 1 : numel(trlinfo)
    trlinfo(j,1) = j;
    trlinfo(j,2) = data_hits.trialinfo{j}.condition;
    trlinfo(j,3) = data_hits.trialinfo{j}.trigger;
    trlinfo(j,4) = data_hits.trialinfo{j}.hit;
end

min_samp = allsub_minsample(s,1);

% now subsample the trials from Cognate and Non-cognate conditions;        
idx_n = find(trlinfo(:,2) == 1); 
idx_nc = find(trlinfo(:,2) == 2); 
rand_trial_c = sortrows(idx_n(randperm(length(idx_n),min_samp)));
rand_trial_nc = sortrows(idx_nc(randperm(length(idx_nc),min_samp)));

% select subsamples trials;
trial_c_temp = data_hits.trial(1,rand_trial_c);
trial_nc_temp = data_hits.trial(1,rand_trial_nc);

% corresponding trial inf;  
info_c_temp = data_hits.trialinfo(rand_trial_c);
info_nc_temp = data_hits.trialinfo(rand_trial_nc); 

% create temporary structure;
data_temp = [];
data_temp.trial = cat(2,trial_c_temp,trial_nc_temp);
data_temp.trialinfo = cat(1,info_c_temp,info_nc_temp);
data_temp.fsample = data_hits.fsample;
data_temp.label = data_hits.label;

% add corresponding time info;
time_temp = data_hits.time{1,1};
for tm = 1: length(data_temp.trialinfo)     
    data_temp.time{1,tm} = time_temp;
end

% calculate Mutual index between theta-Beta at all electrodes over the scalp;
clear pac_beta
cfg = [];
cfg.trials = 'all';
if strcmp(task,'L1N')==1
    cfg.latency = stats_beta.latency_n;
elseif strcmp(task,'L1S')==1
    cfg.latency = stats_beta.latency_s;
end
cfg.phasefreq = [-0.5 0.5] + round(allsubj_meanpeak_info(s,2));
cfg.powerfreq = [-5 5] + round(allsubj_meanpeak_info(s,4));
cfg.nbins = 12;
cfg.keeptrials = 'yes';
cfg.zscore = 'yes';
pac_beta = uni_getBasicPAC(cfg,data_temp);

% save pac data of each participant and the subsampled trials;
fprintf('\nsave pac_gamma in the whole scalp of the subject...\n')
save([pac_out,name2(1:4),name2(5:end-13),'pac_beta'],'pac_beta','-v7.3');
save([pac_out2,name2(1:4),name2(5:end-13),'rand_trials'],'rand_trial_c','rand_trial_nc','-v7.3');
clearvars -except clean_dir tfr_dir tfr_dir pac_out pac_out2 subj_list subj_list_peak task expe ...
                  stats_beta pac_search_light{1, 1}.stats_ttest.p_ttest peak_dir p_value allsub_minsample
                  
end

end

fprintf('\nAll pac_beat in whole brain of the subjects done!\n')

%% Perform Searchlight analysis;

% loop over tasks to compute pac analysis;
expe = {'L1N','L1S'};
close all; pac_anova = []; count2 = 0;
for item = 1 : length(expe)
    
    task = expe{item};
    subj_list2 = dir([pac_dir,'*.mat']); 
    
    % chose electrodes of interest where to average the PAC;
    if strcmp(task,'L1N')==1
        roi_index = find(ismember(all_labels,stats_beta.L1N.roi));
    elseif strcmp(task,'L1S')==1
        roi_index = find(ismember(all_labels,stats_beta.L1S.roi));
    end

    [searchlight_beta,pac_anova] = pac_stats(subj_list2,pac_dir,roi_index,task,pac_anova,count2);
    pac_searchlight{item} = searchlight_beta;
    count2 = count2+2;
    clear searchlight_beta

end

ttest_pac_anova = [];
% T-test mean power against zeros: pvalue(1) t-value(2) and Cohen's d(3);
for ii = 1 : size(pac_anova,2)
    [~,ttest_pac_anova(1,ii),~,stat] = ttest(pac_anova(:,ii),0,'tail','both');
    ttest_pac_anova(2,ii) = stat.tstat;
    ttest_pac_anova(3,ii) = mean(pac_anova(:,ii))/std(pac_anova(:,ii));   
end
% multiple comparison correction of pvalues;
ttest_pac_anova(1,:) = ttest_pac_anova(1,:)*length(ttest_pac_anova(1,:));

% save searchlight pac analysis;
save([pac_dir,'searchlight_theta_beta_pac'], 'pac_searchlight','pac_anova');
close all;

%%
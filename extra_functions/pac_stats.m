function [searchlight_beta,pac_anova] = pac_stats(subj_list,pac_dir,roi_index,task,pac_anova,count2)
    
% select only pac from expe(item);
clear idx
for ii = 1 : length(subj_list)
    idx(1,ii) = contains(subj_list(ii).name,task);
end
idx = find(idx==1);

% cycle through participants;
count = 0;
for s = idx 
  
% define a few more variables and load regressed pac data;
name = subj_list(s).name;
load([pac_dir,name],'pac_beta');
count = count+1;  

% reformat freq
pac_beta = rmfield(pac_beta,{'phase','bin_val','phapow','powpow'});
pac_beta.powspctrm = permute(pac_beta.powspctrm,[4 1 2 3]);
pac_beta.dimord = 'rpt_chan_freq_time';
       
cfg = [];
cfg.keeptrial  = 'no';
cfg.parameter = 'powspctrm';
cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,2), pac_beta.trialinfo, 'UniformOutput', false)));
pac_nc = ft_freqdescriptives(cfg,pac_beta);
cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,1), pac_beta.trialinfo, 'UniformOutput', false)));
pac_c = ft_freqdescriptives(cfg,pac_beta);
    
    % extract PAC data
    for i = 1 : numel(pac_beta)
        group_pac_nc{count,i} = pac_nc;
        group_pac_c{count,i} = pac_c;
    end
    
end

% get grand data structure;
for i = 1
    grand_pac_nc{i} = uni_combineData(group_pac_nc(:,i),{'powspctrm'});
    grand_pac_c{i} = uni_combineData(group_pac_c(:,i),{'powspctrm'});
end

% create dummy freq structure
freq1 = struct('label',{{'pac'}},'time',1,'freq',1,'powspctrm',mean(grand_pac_nc{1}.powspctrm(:,roi_index),2),'dimord','subj_chan','cfg',[]);
freq2 = struct('label',{{'pac'}},'time',1,'freq',1,'powspctrm',mean(grand_pac_c{1}.powspctrm(:,roi_index),2),'dimord','subj_chan','cfg',[]);

% extract metrics
pac_nc = freq1.powspctrm;
pac_c = freq2.powspctrm;

% Dependent t-test to test of PAC between Sync and Async Encoding;
[~,p_ttest,~,stats_ttest] = ttest(pac_nc,pac_c,'tail','both');

% store pac data for anova;
pac_anova(:,[1 2]+ count2) = [pac_nc pac_c];

%% % Searchlight permutation-based approach: compare the original t-value of beta roi against same-size random clusters;

% get pow difference at beta roi; 
freq_orig = freq2;
freq_orig.powspctrm = (freq2.powspctrm - freq1.powspctrm)./freq1.powspctrm;

% set random seed;
rng(1);

% define null hyp (pow diff=zero);
null_hyp = freq_orig;
null_hyp.powspctrm = zeros(size(null_hyp.powspctrm));

% define stat design;
design = zeros(2,size(freq_orig.powspctrm,1)*2);
design(1,:) = repmat(1:size(freq_orig.powspctrm,1),[1 2]);
design(2,:) = [ones(1,size(freq_orig.powspctrm,1)),ones(1,size(freq_orig.powspctrm,1))+1];

% define stat config structure and run ttest to compute t-vlaue of original effect size);
cfg = [];
cfg.method = 'montecarlo';
cfg.correctm = 'no';
cfg.numrandomization = 500;
cfg.ivar = 2;
cfg.uvar = 1;
cfg.parameter = 'powspctrm';
cfg.design = design;
cfg.statistic = 'ft_statfun_depsamplesT';  
cfg.tail = 0;
cfg.correcttail = 'prob';
stat = ft_freqstatistics(cfg,freq_orig,null_hyp);
orig_t = stat;

% put power in a dummy structure;
clear dum_pow
dum_pac = (grand_pac_c{1}.powspctrm - grand_pac_nc{1}.powspctrm)./grand_pac_nc{1}.powspctrm;

% cycle through each random cluster from scalp electrodes (chose the number of permutations with perm = xxx);
tot_perm = 2000; square_t = zeros(tot_perm,1); tot_elec = 1:1:length(grand_pac_c{1,1}.label);

for perm = 1 : tot_perm

    % pick random electrodes to form a cluster with the same size of beta roi;  
    rand_clust = tot_elec(randperm(numel(tot_elec),length(roi_index)));     
    clust_pac = mean(dum_pac(:,rand_clust),2);        
    freq_orig.powspctrm = clust_pac;    

    % get pow;
    stat = ft_freqstatistics(cfg,freq_orig,null_hyp);
    square_t(perm,1) = stat.stat;
    clear rand_clust clust_pow
        
end
    
% get number of instances when ori_t is less than square_t, then get p-value;
count = sum(abs(orig_t.stat) < abs(square_t));
pval = (count+1)/numel(square_t+1);

searchlight_beta.pval = pval;
searchlight_beta.orig_t = orig_t;
searchlight_beta.square_t = square_t;
searchlight_beta.stats_ttest = stats_ttest;
searchlight_beta.stats_ttest.p_ttest = p_ttest;

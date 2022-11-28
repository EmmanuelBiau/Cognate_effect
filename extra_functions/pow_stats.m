function searchlight_beta = pow_stats(L1_nc,L1_c,roi_index,latency,freq,cluster_layout,lay)
    
% reshpae structure to apply data combination;
if size(L1_nc,2) ~= 1
    L1_nc = L1_nc';
    L1_c = L1_c';
end

% get grand data structure;
for i = 1
    grand_pow_c{i} = uni_combineData(L1_nc(:,i),{'powspctrm'});
    grand_pow_nc{i} = uni_combineData(L1_c(:,i),{'powspctrm'});
end

% create dummy freq structures and select time-window/electrodes of interest;
load('layout.mat','lay');
[~,t1]= min(abs(grand_pow_c{1,1}.time-latency(1)));
[~,t2] = min(abs(grand_pow_c{1,1}.time-latency(2)));
freq1 = struct('label',{{'pow'}},'time',1,'freq',1,'powspctrm',mean(mean(mean(grand_pow_nc{1}.powspctrm(:,roi_index,freq(1):freq(2),t1:t2),4),3),2),'dimord','subj_chan','cfg',[]);
freq2 = struct('label',{{'pow'}},'time',1,'freq',1,'powspctrm',mean(mean(mean(grand_pow_c{1}.powspctrm(:,roi_index,freq(1):freq(2),t1:t2),4),3),2),'dimord','subj_chan','cfg',[]);

% extract metrics;
pow_nc = freq1.powspctrm;
pow_c = freq2.powspctrm;
[~,p_ttest,~,stats_ttest] = ttest(pow_nc,pow_c,'tail','both');

%% % Searchlight permutation-based approach: compare the original t-value of beta roi against same-size random clusters;

% get pow difference at beta roi; 
freq_orig = freq1;
freq_orig.powspctrm = (freq1.powspctrm - freq2.powspctrm);

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
clear dum_pow pow
dum_pow = (grand_pow_nc{1}.powspctrm - grand_pow_c{1}.powspctrm);
dum_pow = mean(mean(dum_pow(:,:,freq(1):freq(2),t1:t2),4),3);
clear pow

% cycle through each random cluster from scalp electrodes (chose the number of permutations with perm = xxx);
tot_perm = 2000; square_t = zeros(tot_perm,1); tot_elec = 1:1:length(lay.label)-6;

for perm = 1 : tot_perm
    
    rng('shuffle');
    
    % pick a random electrode and the nearest neighbours to form a mini cluster with the same size of beta roi (elec = 9);  
    rand_elec = randperm(numel(tot_elec),1);

    rand_clust = cluster_layout{rand_elec,5};
        
    if size(rand_clust,1) > size(roi_index,1)
        cut = size(rand_clust,1) - size(roi_index,1);
        rand_clust(end-cut+1:end) = [];
    end

    % get mean pow of the mini cluster;   
    rand_clust = sortrows(rand_clust,'ascend');
    clust_pow = mean(dum_pow(:,rand_clust),2);        
    freq_orig.powspctrm = clust_pow;    

    % perform ttest (0 = two-tailed);
    stat = ft_freqstatistics(cfg,freq_orig,null_hyp);
    square_t(perm,1) = stat.stat;
    clear rand_clust clust_pow rand_x rand_y freq_orig.powspctrm
        
end

% get number of instances when ori_t is less than square_t, then get p-value;
count = sum(abs(square_t) > abs(orig_t.stat));
pval = (count+1)/numel(square_t+1);

searchlight_beta.pval = pval;
searchlight_beta.orig_t = orig_t;
searchlight_beta.square_t = square_t;
searchlight_beta.stats_ttest = stats_ttest;
searchlight_beta.stats_ttest.p_ttest = p_ttest;

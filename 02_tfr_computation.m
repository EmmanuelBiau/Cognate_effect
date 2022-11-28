% compute TFRs of 
clearvars;clc;
clean_dir = 'XXX';
tfr_dir = 'XXX';
tfr_out = 'XXX';

subj_list = dir([clean_dir,'*.mat']);

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
  
    name = subj_list(s).name; load([clean_dir,name]);   

    cfg = [];
    cfg.output = 'pow';
    cfg.method = 'wavelet'; 
    cfg.foi = 1:1:40;
    cfg.toi = -2.5:0.01:2.5;
    cfg.width = 5;
    cfg.keeptrials = 'yes';
    freq_raw = ft_freqanalysis(cfg, data_clean);   
    save([tfr_dir,name(1:3),'_',name(16:end-4),'_rawpow'], 'freq_raw');
    clearvars -except clean_dir subj_list tfr_dir idx_subj expe task

end

end

% Correct power spectrum;
subj_list = dir([tfr_dir,'*.mat']);
data_list = dir([tfr_dir,'*.mat']);

% cycle over tasks;
expe = {'L1N','L1S'};

% select only pac from expe(item);
clear idx
for ii = 1 : length(data_list)
    idx_subj(1,ii) = (contains(data_list(ii).name,expe{1}) || contains(data_list(ii).name,expe{2}));
end
idx_subj = find(idx_subj==1);

for s = idx_subj

    name = subj_list(s).name;
    load([tfr_dir,name]);  

    % correct tfr spectrum;
    freq{1,1} = freq_raw;
    % drop unnecessary fields;
    freq_tmp1 = rmfield(freq{1,1},{'cfg','cumtapcnt'});   
    % subtract 1/f;
    freq_tmp2 = subtract1f(freq_tmp1,freq);
    freq_corr = freq_tmp2;
    clear freq_tmp1 freq_tmp2 freq

    % fprintf('\nsave gamma peak in Hpc of the subject...\n')
    save([tfr_out,name(1:4),name(5:end-11),'_corrpow'],'freq_corr','-v7.3');
    clearvars -except subj_list tfr_dir tfr_out

end

fprintf('\nPower peaks for all subjects done...\n')

% Import tfr data;
subj_list = dir([tfr_dir,'*.mat']);

clear allsubj_L1N_C allsubj_L1N_NC allsubj_L1S_C allsubj_L1S_NC 

% loop over participants;
count_l1n = 0;count_l1s = 0;
for s = 1 : length(subj_list) 

    name = subj_list(s).name;
    load([tfr_dir,name]);   

    if strcmp(freq_corr.trialinfo{1,1}.expe,'L1N')==1  
        count_l1n = count_l1n+1;
        cfg = [];
        cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,1), freq_corr.trialinfo, 'UniformOutput', false)));
        l1c = ft_selectdata(cfg,freq_corr);
        cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,2), freq_corr.trialinfo, 'UniformOutput', false)));
        l1nc = ft_selectdata(cfg,freq_corr);
        allsubj_L1N_C{count_l1n} = l1c;
        allsubj_L1N_NC{count_l1n} = l1nc;
    elseif strcmp(freq_corr.trialinfo{1,1}.expe,'L1S')==1  
        count_l1s = count_l1s+1;
        cfg = [];
        cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,1), freq_corr.trialinfo, 'UniformOutput', false)));
        l1c = ft_selectdata(cfg,freq_corr);
        cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,2), freq_corr.trialinfo, 'UniformOutput', false)));
        l1nc = ft_selectdata(cfg,freq_corr);
        allsubj_L1S_C{count_l1s} = l1c;
        allsubj_L1S_NC{count_l1s} = l1nc;    
    end

end

% save tfr data;
save([tfr_out,'allsubj_tfrs'],'allsubj_L1N_C','allsubj_L1N_NC','allsubj_L1S_C','allsubj_L1S_NC','-v7.3'); 

% average trial within conditions for each participant;
cfg = [];
cfg.keeptrials = 'no';
for ii = 1 : length(allsubj_L1N_C)     
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.hit,1), allsubj_L1N_C{1,ii}.trialinfo, 'UniformOutput', false)));   
    allsubj_L1N_C{ii} = ft_freqdescriptives(cfg,allsubj_L1N_C{ii});
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.hit,1), allsubj_L1N_NC{1,ii}.trialinfo, 'UniformOutput', false)));   
    allsubj_L1N_NC{ii} = ft_freqdescriptives(cfg,allsubj_L1N_NC{ii});
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.hit,1), allsubj_L1S_C{1,ii}.trialinfo, 'UniformOutput', false)));   
    allsubj_L1S_C{ii} = ft_freqdescriptives(cfg,allsubj_L1S_C{ii});
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.hit,1), allsubj_L1S_NC{1,ii}.trialinfo, 'UniformOutput', false)));   
    allsubj_L1S_NC{ii} = ft_freqdescriptives(cfg,allsubj_L1S_NC{ii});     
end

% baseline correction;
cfg = [];
cfg.baseline = [-0.7 -0.2]; 
cfg.baselinetype = 'relchange';
cfg.parameter = 'powspctrm';  
for ii = 1 : length(allsubj_L1N_C)  
    bs_allsubj_L1N_C{ii} = ft_freqbaseline(cfg,allsubj_L1N_C{ii});
    bs_allsubj_L1N_NC{ii} = ft_freqbaseline(cfg,allsubj_L1N_NC{ii});
    bs_allsubj_L1S_C{ii} = ft_freqbaseline(cfg,allsubj_L1S_C{ii});
    bs_allsubj_L1S_NC{ii} = ft_freqbaseline(cfg,allsubj_L1S_NC{ii}); 
end
    
% grand average across participants within conditions;
cfg = [];
cfg.keepindividual = 'yes'; 
tfr_L1N_C = ft_freqgrandaverage(cfg,bs_allsubj_L1N_C{:});
tfr_L1N_NC = ft_freqgrandaverage(cfg,bs_allsubj_L1N_NC{:});
tfr_L1S_C = ft_freqgrandaverage(cfg,bs_allsubj_L1S_C{:});
tfr_L1S_NC = ft_freqgrandaverage(cfg,bs_allsubj_L1S_NC{:});

% Statistic analysis;
load layout;

% prepare structure to store stats analysis;
stats_beta = [];

% prepare design for CBP test;
cfg = [];
nsubj = numel(bs_allsubj_L1N_NC);
design = zeros(2,2*nsubj);
for i = 1:nsubj
  design(1,i) = i;
end
for i = 1:nsubj
  design(1,nsubj+i) = i;
end
design(2,1:nsubj) = 1;
design(2,nsubj+1:2*nsubj) = 2;

% inputs to analyse;
stats_allfreq.latency_allfreq = [-0.1 0.5];
freq = [1 40];
channels = 'all';

% prepare layout for CBP tests;
cfg.method = 'triangulation';
cfg.layout = lay;
nb_elec = ft_prepare_neighbours(cfg);
cfg.channel = channels;
cfg.latency = stats_allfreq.latency_allfreq;
cfg.frequency = freq;
cfg.avgovertime = 'no'; cfg.avgoverchan = 'no'; cfg.avgoverfreq = 'no';
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.tail = 0;
cfg.alpha = 0.05;
cfg.numrandomization = 2000;
cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;
cfg.neighbours = nb_elec;
stats_allfreq.L1N = ft_freqstatistics(cfg,bs_allsubj_L1N_C{:}, bs_allsubj_L1N_NC{:});
stats_allfreq.L1S = ft_freqstatistics(cfg,bs_allsubj_L1S_C{:}, bs_allsubj_L1S_NC{:});
stats_allfreq.tail = cfg.tail;

% prepare design for CBP test;
cfg = [];
nsubj = numel(bs_allsubj_L1N_NC);
design = zeros(2,2*nsubj);
for i = 1:nsubj
  design(1,i) = i;
end
for i = 1:nsubj
  design(1,nsubj+i) = i;
end
design(2,1:nsubj) = 1;
design(2,nsubj+1:2*nsubj) = 2;

% inputs to analyse (beta: 20-35Hz);
stats_beta.latency_n = [0.16 0.26];
stats_beta.latency_s = [0.26 0.38];
stats_beta.freq_beta = [25 35];
channels = 'all';

% prepare layout for CBP tests;
cfg.method = 'triangulation';
cfg.layout = lay;
nb_elec = ft_prepare_neighbours(cfg);
cfg.channel = channels;
cfg.latency = stats_beta.latency_n;
cfg.frequency = stats_beta.freq_beta;
cfg.avgovertime = 'yes'; cfg.avgoverchan = 'no'; cfg.avgoverfreq = 'yes';
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.tail = -1;
cfg.numrandomization = 2000;
cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;
cfg.neighbours = nb_elec;
cfg.minnbchan = 1;
stats_beta.L1N = ft_freqstatistics(cfg,bs_allsubj_L1N_C{:}, bs_allsubj_L1N_NC{:});
stats_beta.L1N.roi = stats_beta.L1N.label(stats_beta.L1N.mask==1);
cfg.latency = stats_beta.latency_s;
stats_beta.L1S = ft_freqstatistics(cfg,bs_allsubj_L1S_C{:}, bs_allsubj_L1S_NC{:});
stats_beta.L1S.roi = stats_beta.L1S.label(stats_beta.L1S.mask==1);
stats_beta.tail = cfg.tail;

% save CBP stats on beta power; 
% save([tfr_dir,'stats_beta'],'stats_beta','stats_allfreq');

% Plot topo of CBP test (check to have the same inputs);
close all; 
cfg = [];
cfg.xlim = stats_beta.latency_n;
cfg.ylim = stats_beta.freq_beta; 
cfg.zlim = [-7 4]; 
cfg.channel = channels;      
cfg.layout = lay; 
cfg.marker = 'off';
cfg.gridscale = 300;                  
cfg.style = 'fill';  
cfg.parameter = 'stat';
cfg.mask = 'mask';
cfg.highlight = 'on';
cfg.comment = 'no';
cfg.highlightchannel = stats_beta.L1N.label(stats_beta.L1N.mask==1);
s1 = subplot(1,2,1); ft_topoplotTFR(cfg,stats_beta.L1N);hold on
caxis([-4 4]); axis xy; title('L1 Naming Non-Cognate > Cognate','position',[0, 0.8, 0]);
c1 = colorbar(s1,'manual','southoutside','position',[0.1681 0.1740 0.2581 0.0508],'ylim', [-4 4]);
cfg.xlim = stats_beta.latency_s;
cfg.highlightchannel = stats_beta.L1S.label(stats_beta.L1S.mask==1);
s2 = subplot(1,2,2); ft_topoplotTFR(cfg,stats_beta.L1S);
caxis([-4 4]); axis xy; title('L1 Semantic Non-Cognate > Cognate','position',[0, 0.8, 0]);
c2 = colorbar(s2,'southoutside','position',[0.61 0.1740 0.2581 0.0508],'ylim', [-4 4]);
set(gcf,'color','w','position',[266.6 226.6 904 420]);

% Plot TFR of T-values all frequencies at pool of significant electrodes;
cfg = [];
cfg.xlim = [-0.1 0.5];
cfg.ylim = [1 40];  
cfg.zlim = [-5 2]; 
cfg.layout = lay;
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskalpha = 0.8;
cfg.colorbar ='yes';
cfg.style = 'fill';  

figure; 
subplot(1,2,1);
cfg.channel = stats_beta.L1N.label(stats_beta.L1N.mask==1);
ft_singleplotTFR(cfg,stats_allfreq.L1N);
set(gcf,'color','w');
caxis([-3 3]); axis xy; xlabel('Time [s]');ylabel('Frequencies [Hz]');
title('L1 Naming Non-Cognate > Cognate');
subplot(1,2,2);
cfg.channel = stats_beta.L1S.label(stats_beta.L1S.mask==1);
ft_singleplotTFR(cfg,stats_allfreq.L1S);
caxis([-3 3]); axis xy; xlabel('Time [s]');ylabel('Frequencies [Hz]');
title('L1 Semantic Non-Cognate > Cognate');
set(gcf,'color','w','position',[266.6 226.6 904 420]);

%% Perform t-test and searchlight analysis on beta power;

tfr_searchlight = [];
tfr_searchlight.L1N_C = bs_allsubj_L1N_C;
tfr_searchlight.L1N_NC = bs_allsubj_L1N_NC;
tfr_searchlight.L1S_C = bs_allsubj_L1S_C;
tfr_searchlight.L1S_NC = bs_allsubj_L1S_NC;

% prepare raw mini_clusters;
cfg = [];
cfg.method = 'distance';
cfg.layout = lay;
mini_clusts = ft_prepare_miniclust(cfg);

cluster_layout = [];
for ii = 1: size(mini_clusts,2)
    cluster_layout{ii,1} = mini_clusts(ii).label;
    cluster_layout{ii,2} = mini_clusts(ii).distance;
    cluster_layout{ii,3} = mini_clusts(ii).neighblabel;
end

idx_remv = [];
idx_remv = find(strcmp(cluster_layout(:,1),'EOG1')==1 | strcmp(cluster_layout(:,1),'EOG2')==1 | strcmp(cluster_layout(:,1),'MSD1')==1 | strcmp(cluster_layout(:,1),'MSD2')==1);
cluster_layout(idx_remv,:) = []; clear idx_remv

% sort electrodes by their relative distance to the center of the mini cluster;
for ii = 1 : size(cluster_layout,1)
    
    [cluster_layout{ii,2},idx_dist] = sortrows(cluster_layout{ii,2},'ascend');
    cluster_layout{ii,3} = cluster_layout{ii,3}(idx_dist);

    idx_remv = [];
    idx_remv = find(strcmp(cluster_layout{ii,3},'EOG1')==1 | strcmp(cluster_layout{ii,3},'EOG2')==1 | strcmp(cluster_layout{ii,3},'MSD1')==1 | strcmp(cluster_layout{ii,3},'MSD2')==1);
    cluster_layout{ii,2}(idx_remv) = []; cluster_layout{ii,3}(idx_remv) = [];  
    clear idx_remv

    cluster_layout{ii,4} = [cluster_layout{ii,1};cluster_layout{ii,3}];
    [~,cluster_layout{ii,5}] = ismember(cluster_layout{ii,4},lay.label);
    clear cut 
    
end
clc;

% loop over tasks to compute pac analysis;
searchlight_beta_pow = {};
expe = {'L1N','L1S'}; searchlight_beta_pow = {};
for item = 1 : length(expe)
    
    task = expe{item};

    % chose frequency band of interest;
    freq = [25 35];

    % chose electrodes/latency of interest where to average the power;
    if strcmp(task,'L1N')==1   
        L1_nc = tfr_searchlight.L1N_NC; L1_c = tfr_searchlight.L1N_C; 
        roi_index = find(ismember(all_labels,stats_beta.L1N.roi));
        latency = stats_beta.latency_n;
    elseif strcmp(task,'L1S')==1
        L1_nc = tfr_searchlight.L1S_NC; L1_c = tfr_searchlight.L1S_C;
        roi_index = find(ismember(all_labels,stats_beta.L1S.roi));
        latency = stats_beta.latency_s;
    end

    searchlight_beta = pow_stats(L1_nc,L1_c,roi_index,latency,freq,cluster_layout,lay);
    searchlight_beta_pow{item} = searchlight_beta;
    clear searchlight_beta

end

save([tfr_dir,'searchlight_beta_pow'],'searchlight_beta_pow');

% plot the distribution of t-values of the 2000 permutations and the original data;
close all; figure; 
subplot(1,2,1);
h = histogram(searchlight_beta_pow{1,1}.square_t,10,'Normalization','probability');hold on
plot([round(searchlight_beta_pow{1,1}.orig_t.stat,2) round(searchlight_beta_pow{1,1}.orig_t.stat,2)],[0 0.25],':k','Linewidth',2);
xlabel('t-values');ylabel('Searchlight prob.');
ylim([0 0.25]); xlim([-4 1]); xticks([-4 -2 0]);
subplot(1,2,2);
h = histogram(searchlight_beta_pow{1,2}.square_t,8,'Normalization','probability');hold on
plot([round(searchlight_beta_pow{1,2}.orig_t.stat,2) round(searchlight_beta_pow{1,2}.orig_t.stat,2)],[0 0.25],':k','Linewidth',2);
xlabel('t-values');ylabel('Searchlight prob.');
ylim([0 0.25]); xlim([-5 -1]); xticks([-5 -3 -1]);
set(gcf,'position',[366.6 204.2 647.2 352.8]);

%%
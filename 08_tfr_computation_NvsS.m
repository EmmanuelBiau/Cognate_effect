%% Put everyone with together;
clearvars;clc; 
tfr_dir = 'XXX'; load([tfr_dir,'allsubj_tfrs']);

clear allsubj_L1N allsubj_L1S
for ii = 1 : length(allsubj_L1N_NC)
    allsubj_L1N{ii} = allsubj_L1N_NC{ii};
    allsubj_L1N{ii}.powspctrm = cat(1,allsubj_L1N_NC{ii}.powspctrm,allsubj_L1N_C{ii}.powspctrm);
    allsubj_L1N{ii}.trialinfo = [allsubj_L1N{ii}.trialinfo;allsubj_L1N_C{ii}.trialinfo];
    allsubj_L1S{ii} = allsubj_L1S_NC{ii};
    allsubj_L1S{ii}.powspctrm = cat(1,allsubj_L1S_NC{ii}.powspctrm,allsubj_L1S_C{ii}.powspctrm);
    allsubj_L1S{ii}.trialinfo = [allsubj_L1S{ii}.trialinfo;allsubj_L1S_C{ii}.trialinfo];
end

% average trial within conditions for each participant;
cfg = [];
cfg.keeptrials = 'no';
for ii = 1 : length(allsubj_L1N)     
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.hit,1), allsubj_L1N{1,ii}.trialinfo, 'UniformOutput', false)));   
    allsubj_L1N{ii} = ft_freqdescriptives(cfg,allsubj_L1N{ii});
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.hit,1), allsubj_L1S{1,ii}.trialinfo, 'UniformOutput', false)));   
    allsubj_L1S{ii} = ft_freqdescriptives(cfg,allsubj_L1S{ii});
end

% baseline correction;
cfg = [];
cfg.baseline = [-0.7 -0.2]; 
cfg.baselinetype = 'relchange';
cfg.parameter = 'powspctrm';  
for ii = 1 : length(allsubj_L1N)  
    bs_allsubj_L1N{ii} = ft_freqbaseline(cfg,allsubj_L1N{ii});
    bs_allsubj_L1S{ii} = ft_freqbaseline(cfg,allsubj_L1S{ii});
end

% grand average across participants within conditions;
cfg = [];
cfg.keepindividual = 'yes'; 
tfr_L1N_C = ft_freqgrandaverage(cfg,bs_allsubj_L1N{:});
tfr_L1N_NC = ft_freqgrandaverage(cfg,bs_allsubj_L1S{:});

% Compute TFRs with t-tvalues;
clc; load layout;

% prepare design for CBP test;
cfg = [];
nsubj = numel(allsubj_L1N);
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
latency_allfreq = [-0.1 0.5];
freq = [1 40];
channels = 'all';

% prepare layout for CBP tests;
clear stats_allfreq
cfg.method = 'triangulation';
cfg.layout = lay;
nb_elec = ft_prepare_neighbours(cfg);
cfg.channel = channels;
cfg.latency = latency_allfreq;
cfg.frequency = freq;
cfg.avgovertime = 'no'; cfg.avgoverchan = 'no'; cfg.avgoverfreq = 'no';
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 2000;
cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;
cfg.neighbours = nb_elec;
stats_allfreq_NvsS = ft_freqstatistics(cfg,bs_allsubj_L1N{:},bs_allsubj_L1S{:});
stats_allfreq_NvsS.tail = cfg.tail;
stats_allfreq_NvsS.latency_allfreq = [-0.1 0.5];

% save stats all frequencies; 
% save([tfr_dir,'stats_NvsS'],'stats_allfreq_NvsS');

% Plot TFR of T-values all frequencies at pool of significant electrodes;
elec_com = {'C3','CP3','P3','C4','CP4','P4','Cz','CPz','Pz'}; 

close all;
cfg = [];
cfg.xlim = [-0.1 0.5];
cfg.ylim = [2 40];  
cfg.layout = lay;
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskalpha = 0;
cfg.colorbar ='yes';
cfg.style = 'fill';  
close all; figure; 
cfg.channel = elec_com;
ft_singleplotTFR(cfg,stats_allfreq_NvsS);
set(gcf,'color','w');
caxis([-4 4]); axis xy; xlabel('Time [s]');ylabel('Frequencies [Hz]');
title('L1 Naming vs Semantic');

%% 
%% Put everyone with together;
clearvars;clc;
tfr_dir = 'XXX'; load([tfr_dir,'allsubj_tfrs']); load layout;

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

stats_theta_L1NvsS = [];

% prepare design for CBP test;
cfg = [];
nsubj = numel(bs_allsubj_L1N_C);
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
stats_theta_L1NvsS.latency_n = [0.16 0.26];
stats_theta_L1NvsS.latency_s = [0.26 0.38];
stats_theta_L1NvsS.freq_beta = [3 7];
channels = 'all';

% prepare layout for CBP tests;
cfg.method = 'triangulation';
cfg.layout = lay;
nb_elec = ft_prepare_neighbours(cfg);
cfg.channel = channels;
cfg.latency = stats_theta_L1NvsS.latency_n;
cfg.frequency = stats_theta_L1NvsS.freq_beta;
cfg.avgovertime = 'yes'; cfg.avgoverchan = 'no'; cfg.avgoverfreq = 'yes';
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.tail = 0;
cfg.numrandomization = 2000;
cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;
cfg.neighbours = nb_elec;
cfg.minnbchan = -11;
stats_theta_L1NvsS = ft_freqstatistics(cfg,bs_allsubj_L1N_C{:}, bs_allsubj_L1S_C{:});
stats_theta_L1NvsS.roi = stats_theta_L1NvsS.label(stats_theta_L1NvsS.mask==1);
stats_theta_L1NvsS.tail = cfg.tail;

close all; 
cfg = [];
cfg.channel = 'all';      
cfg.layout = lay; 
cfg.marker = 'off';
cfg.gridscale = 300;                  
cfg.style = 'fill';  
cfg.parameter = 'stat';
cfg.comment = 'no';
s1 = figure; ft_topoplotTFR(cfg,stats_theta_L1NvsS);
caxis([-4 4]); axis xy; title('Size-judgment > Picture-naming','position',[0,0.8,0]);
colorbar; caxis([-4 4]);
set(gcf,'color','w');

%%
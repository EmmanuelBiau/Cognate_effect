clearvars;clc;
tfr_dir = 'XXX'; 
clean_dir = 'XXX';
erp_dir = 'XXX';

load([tfr_dir,'stats_beta.mat']); load layout;

subj_list = dir([clean_dir,'*.mat']);

% loop over tasks;
expe = {'L1N','L1S'};
close all;

% loop over participants;
count_l1n = 0;count_l1s = 0;
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
    cfg.latency = [-0.2 1];
    data_clean = ft_selectdata(cfg,data_clean);

    if strcmp(data_clean.trialinfo{1,1}.expe,'L1N')==1  
        count_l1n = count_l1n+1;
        cfg = [];
        cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,1), data_clean.trialinfo, 'UniformOutput', false)));
        l1c = ft_selectdata(cfg,data_clean);
        cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,2), data_clean.trialinfo, 'UniformOutput', false)));
        l1nc = ft_selectdata(cfg,data_clean);
        allsubj_data.L1N_C{count_l1n} = l1c;
        allsubj_data.L1N_NC{count_l1n} = l1nc;
    elseif strcmp(data_clean.trialinfo{1,1}.expe,'L1S')==1  
        count_l1s = count_l1s+1;
        cfg = [];
        cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,1), data_clean.trialinfo, 'UniformOutput', false)));
        l1c = ft_selectdata(cfg,data_clean);
        cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,2), data_clean.trialinfo, 'UniformOutput', false)));
        l1nc = ft_selectdata(cfg,data_clean);
        allsubj_data.L1S_C{count_l1s} = l1c;
        allsubj_data.L1S_NC{count_l1s} = l1nc;    
    end       

end

end

% average trial within conditions for each participant;
cfg = [];
cfg.keeptrials = 'no';
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.2 0];

for ii = 1 : length(allsubj_data.L1N_C)     
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.hit,1), allsubj_data.L1N_C{1,ii}.trialinfo, 'UniformOutput', false)));   
    allsubj_L1N_C{ii} = ft_timelockanalysis(cfg,allsubj_data.L1N_C{ii});
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.hit,1), allsubj_data.L1N_NC{1,ii}.trialinfo, 'UniformOutput', false)));   
    allsubj_L1N_NC{ii} = ft_timelockanalysis(cfg,allsubj_data.L1N_NC{ii});
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.hit,1), allsubj_data.L1S_C{1,ii}.trialinfo, 'UniformOutput', false)));   
    allsubj_L1S_C{ii} = ft_timelockanalysis(cfg,allsubj_data.L1S_C{ii});
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.hit,1), allsubj_data.L1S_NC{1,ii}.trialinfo, 'UniformOutput', false)));   
    allsubj_L1S_NC{ii} = ft_timelockanalysis(cfg,allsubj_data.L1S_NC{ii});     
end

%% Statistic analysis on ERPs over the scalp;
stats_allerp = [];

% prepare design for CBP test;
cfg = [];
nsubj = numel(allsubj_L1N_C);
design = zeros(2,2*nsubj);
for i = 1:nsubj
  design(1,i) = i;
end
for i = 1:nsubj
  design(1,nsubj+i) = i;
end
design(2,1:nsubj) = 1;
design(2,nsubj+1:2*nsubj) = 2;

% prepare layout for CBP tests;
cfg.method = 'triangulation';
cfg.layout = lay;
nb_elec = ft_prepare_neighbours(cfg);
cfg.channel = 'all';
cfg.latency = [0.16 0.26];
cfg.avgovertime = 'yes'; cfg.avgoverchan = 'no';
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.025;
cfg.tail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 2000;
cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;
cfg.neighbours = nb_elec;
stats_allerp.L1N_earlywind = ft_timelockstatistics(cfg,allsubj_L1N_C{:},allsubj_L1N_NC{:});
stats_allerp.L1S_earlywind = ft_timelockstatistics(cfg,allsubj_L1S_C{:},allsubj_L1S_NC{:});
cfg.latency = [0.26 0.38];
stats_allerp.L1N_latewind = ft_timelockstatistics(cfg,allsubj_L1N_C{:},allsubj_L1N_NC{:});
stats_allerp.L1S_latewind = ft_timelockstatistics(cfg,allsubj_L1S_C{:},allsubj_L1S_NC{:});

% save([tfr_dir,'stats_erps.mat'],'stats_allerp');

close all; figure;
cfg = [];
cfg.layout = lay;
cfg.marker = 'off';
cfg.gridscale = 300;                  
cfg.style = 'fill';  
cfg.parameter = 'stat';
cfg.mask = 'mask';
cfg.comment = 'no';
cfg.highlight = 'on';
cfg.highlightchannel = stats_allerp.L1N_earlywind.label(stats_allerp.L1N_earlywind.mask==1);
s1 = subplot(1,2,1); ft_topoplotER(cfg,stats_allerp.L1N_earlywind);
caxis([-4 4]); axis xy; title('L1 Naming Non-Cognate > Cognate','position',[0, 0.8, 0]);
c1 = colorbar(s1,'manual','southoutside','position',[0.1681 0.1740 0.2581 0.0508],'ylim', [-4 4]);
hold on
cfg.highlightchannel = stats_allerp.L1S_earlywind.label(stats_allerp.L1S_earlywind.mask==1);
s2= subplot(1,2,2); ft_topoplotER(cfg,stats_allerp.L1S_earlywind);
caxis([-4 4]); axis xy; title('L1 Semantic Non-Cognate > Cognate','position',[0, 0.8, 0]);
c2 = colorbar(s2,'southoutside','position',[0.61 0.1740 0.2581 0.0508],'ylim', [-4 4]);
hold on

figure;
cfg = [];
cfg.layout = lay;
cfg.marker = 'off';
cfg.gridscale = 300;                  
cfg.style = 'fill';  
cfg.parameter = 'stat';
cfg.mask = 'mask';
cfg.comment = 'no';
cfg.highlight = 'on';
cfg.highlightchannel = stats_allerp.L1N_latewind.label(stats_allerp.L1N_latewind.mask==1);
s1 = subplot(1,2,1); ft_topoplotER(cfg,stats_allerp.L1N_latewind);
caxis([-4 4]); axis xy; title('L1 Naming Non-Cognate > Cognate','position',[0, 0.8, 0]);
c1 = colorbar(s1,'manual','southoutside','position',[0.1681 0.1740 0.2581 0.0508],'ylim', [-4 4]);
hold on
cfg.highlightchannel = stats_allerp.L1S_latewind.label(stats_allerp.L1S_latewind.mask==1);
s2= subplot(1,2,2); ft_topoplotER(cfg,stats_allerp.L1S_latewind);
caxis([-4 4]); axis xy; title('L1 Semantic Non-Cognate > Cognate','position',[0, 0.8, 0]);
c2 = colorbar(s2,'southoutside','position',[0.61 0.1740 0.2581 0.0508],'ylim', [-4 4]);

%% Plot ERPs at the rois;

% compute grand averaged ERPs in the four conditions;
cfg = [];
cfg.channel = 'all';
cfg.latency = 'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'no'; 
GA_L1N_C = ft_timelockgrandaverage(cfg,allsubj_L1N_C{:});
GA_L1N_NC = ft_timelockgrandaverage(cfg,allsubj_L1N_NC{:});
GA_L1S_C = ft_timelockgrandaverage(cfg,allsubj_L1S_C{:});
GA_L1S_NC = ft_timelockgrandaverage(cfg,allsubj_L1S_NC{:});

cfg = [];
cfg.layout = lay;
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfreq = 20;
cfg.xlim = [-0.1 0.5];
cfg.channel = stats_allerp.L1N_earlywind.label(stats_allerp.L1N_earlywind.mask==1);
cfg.title = 'Picture-naming';
cfg.linewidth = 2;
close all; figure; subplot(1,2,1);
ft_singleplotER(cfg,GA_L1N_C,GA_L1N_NC);
xticks([-0.1 0 0.1 0.2 0.3 0.4 0.5]); xlim([-0.1 0.5]); %ylim([-1 3]);
hold on; plot([-0.1 0.5],[0 0],':k');
xlabel('Time [s]'); ylabel('Amplitude [uV]');
legend({'Cognate', 'Non-cognate'},'Location','NorthWest')
hold on
ft_singleplotER(cfg,GA_L1S_C,GA_L1S_NC);
xticks([-0.1 0 0.1 0.2 0.3 0.4 0.5]); xlim([-0.1 0.5]); ylim([-2 8]);
hold on; plot([-0.1 0.5],[0 0],':k');
xlabel('Time [s]'); ylabel('Amplitude [uV]');
legend({'Cognate', 'Non-cognate'},'Location','NorthWest')
set(gcf,'Position',[303.4 260.2 682.4 420],'color','w');

%%
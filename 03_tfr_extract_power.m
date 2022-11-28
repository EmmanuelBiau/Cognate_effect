%% Extract the amplitudes for each participants for your time-window/electrodes of interest;
clearvars;clc;
tfr_dir = 'XXX'; load([tfr_dir,'allsubj_tfrs']);load([tfr_dir,'stats_beta']);
behav_out = 'XXX'; load([behav_out,'allsubj_behav']);

% average trial within conditions for each participant;
for ii = 1 : length(allsubj_L1N_C)  
    cfg = [];
    cfg.keeptrials = 'no';
    allsubj_L1N_C{ii} = ft_freqdescriptives(cfg,allsubj_L1N_C{ii});
    allsubj_L1N_NC{ii} = ft_freqdescriptives(cfg,allsubj_L1N_NC{ii});
    allsubj_L1S_C{ii} = ft_freqdescriptives(cfg,allsubj_L1S_C{ii});
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

% parameters;
% frequency ranges;
foi_theta = [3 7];
foi_alpha = [8 12];
foi_beta = [25 35];

% time-windows;
toi_n = stats_beta.latency_n;
toi_s = stats_beta.latency_s;

% rois;
chan_n = stats_beta.L1N.roi;
chan_s = stats_beta.L1S.roi;

% extract beta power;
cfg = [];
cfg.latency = toi_n; cfg.frequency = foi_theta; cfg.channel = chan_n;
cfg.avgovertime = 'yes'; cfg.avgoverfreq = 'yes'; cfg.avgoverchan = 'yes';
theta_pow.L1N_C = ft_selectdata(cfg,tfr_L1N_C); theta_pow.L1N_NC = ft_selectdata(cfg,tfr_L1N_NC);
theta_pow.L1N_C = theta_pow.L1N_C.powspctrm; theta_pow.L1N_NC = theta_pow.L1N_NC.powspctrm;
cfg = [];
cfg.latency = toi_s; cfg.frequency = foi_theta; cfg.channel = chan_s; 
cfg.avgovertime = 'yes'; cfg.avgoverfreq = 'yes'; cfg.avgoverchan = 'yes';
theta_pow.L1S_C = ft_selectdata(cfg,tfr_L1S_C); theta_pow.L1S_NC = ft_selectdata(cfg,tfr_L1S_NC);
theta_pow.L1S_C = theta_pow.L1S_C.powspctrm; theta_pow.L1S_NC = theta_pow.L1S_NC.powspctrm;
% extract alpha power;
cfg = [];
cfg.latency = toi_n; cfg.frequency = foi_alpha; cfg.channel = chan_n;
cfg.avgovertime = 'yes'; cfg.avgoverfreq = 'yes'; cfg.avgoverchan = 'yes';
alpha_pow.L1N_C = ft_selectdata(cfg,tfr_L1N_C); alpha_pow.L1N_NC = ft_selectdata(cfg,tfr_L1N_NC);
alpha_pow.L1N_C = alpha_pow.L1N_C.powspctrm; alpha_pow.L1N_NC = alpha_pow.L1N_NC.powspctrm;
cfg = [];
cfg.latency = toi_s; cfg.frequency = foi_alpha; cfg.channel = chan_s; 
cfg.avgovertime = 'yes'; cfg.avgoverfreq = 'yes'; cfg.avgoverchan = 'yes';
alpha_pow.L1S_C = ft_selectdata(cfg,tfr_L1S_C); alpha_pow.L1S_NC = ft_selectdata(cfg,tfr_L1S_NC);
alpha_pow.L1S_C = alpha_pow.L1S_C.powspctrm; alpha_pow.L1S_NC = alpha_pow.L1S_NC.powspctrm;
% extract beta power;
cfg = [];
cfg.latency = toi_n; cfg.frequency = foi_beta; cfg.channel = chan_n;
cfg.avgovertime = 'yes'; cfg.avgoverfreq = 'yes'; cfg.avgoverchan = 'yes';
beta_pow.L1N_C = ft_selectdata(cfg,tfr_L1N_C); beta_pow.L1N_NC = ft_selectdata(cfg,tfr_L1N_NC);
beta_pow.L1N_C = beta_pow.L1N_C.powspctrm; beta_pow.L1N_NC = beta_pow.L1N_NC.powspctrm;
cfg = [];
cfg.latency = toi_s; cfg.frequency = foi_beta; cfg.channel = chan_s;
cfg.avgovertime = 'yes'; cfg.avgoverfreq = 'yes'; cfg.avgoverchan = 'yes';
beta_pow.L1S_C = ft_selectdata(cfg,tfr_L1S_C); beta_pow.L1S_NC = ft_selectdata(cfg,tfr_L1S_NC);
beta_pow.L1S_C = beta_pow.L1S_C.powspctrm; beta_pow.L1S_NC = beta_pow.L1S_NC.powspctrm;

% T-tests mean power against zeros: pvalue(1) t-value(2) and Cohen's d(3);
clear pow ttest_pow
pow.theta(:,1) = theta_pow.L1N_NC; pow.theta(:,2) = theta_pow.L1N_C; pow.theta(:,3) = theta_pow.L1S_NC; pow.theta(:,4) = theta_pow.L1S_C;
pow.alpha(:,1) = alpha_pow.L1N_NC; pow.alpha(:,2) = alpha_pow.L1N_C; pow.alpha(:,3) = alpha_pow.L1S_NC; pow.alpha(:,4) = alpha_pow.L1S_C;
pow.beta(:,1) = beta_pow.L1N_NC; pow.beta(:,2) = beta_pow.L1N_C; pow.beta(:,3) = beta_pow.L1S_NC; pow.beta(:,4) = beta_pow.L1S_C;

for ii = 1 : size(pow.theta,2)
    
    [~,ttest_pow.theta(1,ii),~,stat] = ttest(pow.theta(:,ii),0,'tail','both');
    ttest_pow.theta(2,ii) = stat.tstat;
    ttest_pow.theta(3,ii) = mean(pow.theta(:,ii))/std(pow.theta(:,ii));
    ttest_pow.theta(4,ii) = mean(pow.theta(:,ii));
    ttest_pow.theta(5,ii) = std(pow.theta(:,ii));
    [~,ttest_pow.alpha(1,ii),~,stat] = ttest(pow.alpha(:,ii),0,'tail','both'); 
    ttest_pow.alpha(2,ii) = stat.tstat;
    ttest_pow.alpha(3,ii) = mean(pow.alpha(:,ii))/std(pow.alpha(:,ii));
    ttest_pow.alpha(4,ii) = mean(pow.alpha(:,ii));
    ttest_pow.alpha(5,ii) = std(pow.alpha(:,ii));
    [~,ttest_pow.beta(1,ii),~,stat] = ttest(pow.beta(:,ii),0,'tail','both');
    ttest_pow.beta(2,ii) = stat.tstat;
    ttest_pow.beta(3,ii) = mean(pow.beta(:,ii))/std(pow.beta(:,ii));
    ttest_pow.beta(4,ii) = mean(pow.beta(:,ii));
    ttest_pow.beta(5,ii) = std(pow.beta(:,ii));   

end

% multiple comparison correction of p-values;
ttest_pow.theta(1,:) = ttest_pow.theta(1,:)*length(ttest_pow.theta(1,:));
ttest_pow.alpha(1,:) = ttest_pow.alpha(1,:)*length(ttest_pow.alpha(1,:));
ttest_pow.beta(1,:) = ttest_pow.beta(1,:)*length(ttest_pow.beta(1,:));

% Correlations between Behavior and Beta power;
clear corr_scores 
corr_scores(:,1) = 1:length(bs_allsubj_L1N_C);
corr_scores(:,2) = normalize(cell2mat(allsubj.nc_cr(strcmp(allsubj.task,'L1N')==1)) - cell2mat(allsubj.c_cr(strcmp(allsubj.task,'L1N')==1)),'zscore');
corr_scores(:,3) = normalize(cell2mat(allsubj.nc_cr(strcmp(allsubj.task,'L1S')==1)) - cell2mat(allsubj.c_cr(strcmp(allsubj.task,'L1S')==1)),'zscore');
corr_scores(:,4) = normalize(cell2mat(allsubj.nc_rt(strcmp(allsubj.task,'L1N')==1)) - cell2mat(allsubj.c_rt(strcmp(allsubj.task,'L1N')==1)),'zscore');
corr_scores(:,5) = normalize(cell2mat(allsubj.nc_rt(strcmp(allsubj.task,'L1S')==1)) - cell2mat(allsubj.c_rt(strcmp(allsubj.task,'L1S')==1)),'zscore');
corr_scores(:,6) = normalize(beta_pow.L1N_NC - beta_pow.L1N_C,'zscore');
corr_scores(:,7) = normalize(beta_pow.L1S_NC - beta_pow.L1S_C,'zscore');
corr_scores(:,8) = normalize(theta_pow.L1N_NC - theta_pow.L1N_C,'zscore');
corr_scores(:,9) = normalize(theta_pow.L1S_NC - theta_pow.L1S_C,'zscore');
corr_scores(:,10) = normalize(alpha_pow.L1N_NC - alpha_pow.L1N_C,'zscore');
corr_scores(:,11) = normalize(alpha_pow.L1S_NC - alpha_pow.L1S_C,'zscore');

%% Plot correlations between behavioural performances and power;
fig = figure;
subplot(1,2,1);
s = scatter(corr_scores(:,6),corr_scores(:,2));
fit = lsline;
s.MarkerFaceColor = [0.5 0.5 0.5];
s.MarkerEdgeColor = [0 0 0];
s.MarkerFaceAlpha = 1;
s.MarkerEdgeAlpha = 1;
s.SizeData = 15;
fit.Color = [0 0 0];
fit.LineWidth = 2;
xlim([-3 3]); ylim([-2 2]);
set(gca,'xtick',[-2 0 2],'ytick',[-2 -1 0 1 2])
[roh(1,1),pval(1,1)] = corr(corr_scores(:,6),corr_scores(:,2),'tail','both','type','Pearson');
text(1,-1.3,['coef = ',num2str(round(roh(1,1),3))],'color','k');
text(1,-1.5,['pval = ',num2str(round(pval(1,1),3))],'color','k');
xlabel('Diff. Beta power (z-scores)'); ylabel('Diff. Accuracy (z-scores)');
title('L1 Naming');
subplot(1,2,2);
s = scatter(corr_scores(:,7),corr_scores(:,3));
fit = lsline;
s.MarkerFaceColor = [0.5 0.5 0.5];
s.MarkerEdgeColor = [0 0 0];
s.MarkerFaceAlpha = 1;
s.MarkerEdgeAlpha = 1;
s.SizeData = 15;
fit.Color = [0 0 0];
fit.LineWidth = 2;
xlim([-3 3]); ylim([-2 2]);
set(gca,'xtick',[-2 0 2],'ytick',[-2 -1 0 1 2])
[roh(1,2),pval(1,2)] = corr(corr_scores(:,7),corr_scores(:,3),'tail','both','type','Pearson');
text(-2.5,-1.3,['coef = ',num2str(round(roh(1,2),3))],'color','k');
text(-2.5,-1.5,['pval = ',num2str(round(pval(1,2),3))],'color','k');
xlabel('Diff. Beta power (z-scores)'); ylabel('Diff. Accuracy (z-scores)');
title('L1 Semantic');
set(fig,'position',[488.2 341.8 711.2 420],'color','w');

[Fisher(1,1),~,~,~] = corr_rtest(roh(1,1),roh(1,2),18,18);

fig2 = figure;
subplot(1,2,1);
s = scatter(corr_scores(:,4),corr_scores(:,6));
fit = lsline;
s.MarkerFaceColor = [0.5 0.5 0.5];
s.MarkerEdgeColor = [0 0 0];
s.MarkerFaceAlpha = 1;
s.MarkerEdgeAlpha = 1;
s.SizeData = 15;
fit.Color = [0 0 0];
fit.LineWidth = 2;
xlim([-3 5]); ylim([-2 2]);
set(gca,'xtick',[-2 0 2 4],'ytick',[-2 -1 0 1 2])
[roh(1,1),pval(1,1)] = corr(corr_scores(:,4),corr_scores(:,6),'tail','both');
text(2,-1.3,['coef = ',num2str(round(roh(1,1),3))],'color','k');
text(2,-1.5,['pval = ',num2str(round(pval(1,1),3))],'color','k');
xlabel('Diff. Reaction Times (z-scores)'); ylabel('Diff. Beta power (z-scores)');
title('L1 Naming');
subplot(1,2,2);
s = scatter(corr_scores(:,5),corr_scores(:,7));
fit = lsline;
s.MarkerFaceColor = [0.5 0.5 0.5];
s.MarkerEdgeColor = [0 0 0];
s.MarkerFaceAlpha = 1;
s.MarkerEdgeAlpha = 1;
s.SizeData = 15;
fit.Color = [0 0 0];
fit.LineWidth = 2;
xlim([-3 5]); ylim([-2 2]);
set(gca,'xtick',[-2 0 2 4],'ytick',[-2 -1 0 1 2])
[roh(1,2),pval(1,2)] = corr(corr_scores(:,5),corr_scores(:,7),'tail','right');
text(2,-1.3,['coef = ',num2str(round(roh(1,2),3))],'color','k');
text(2,-1.5,['pval = ',num2str(round(pval(1,2),3))],'color','k');
xlabel('Diff. Reaction Times (z-scores)'); ylabel('Diff. Beta power (z-scores)');
title('L1 Semantic');
set(fig2,'position',[488.2 341.8 711.2 420],'color','w');

[Fisher(1,2),~,~,~] = corr_rtest(roh(1,1),roh(1,2),18,18);

fig3 = figure;
subplot(1,2,1);
s = scatter(corr_scores(:,2),corr_scores(:,3));
fit = lsline;
s.MarkerFaceColor = [0.5 0.5 0.5];
s.MarkerEdgeColor = [0 0 0];
s.MarkerFaceAlpha = 1;
s.MarkerEdgeAlpha = 1;
s.SizeData = 15;
fit.Color = [0 0 0];
fit.LineWidth = 2;
xlim([-3 5]); ylim([-2 2]);
set(gca,'xtick',[-2 0 2 4],'ytick',[-2 -1 0 1 2])
[roh(1,1),pval(1,1)] = corr(corr_scores(:,2),corr_scores(:,3),'tail','right');
text(1,-1.3,['coef = ',num2str(round(roh(1,1),3))],'color','k');
text(1,-1.5,['pval = ',num2str(round(pval(1,1),3))],'color','k');
xlabel('NC-C Naming (z-scores)'); ylabel('NC-C Semantic (z-scores)');
title('Accuracy');
subplot(1,2,2);
s = scatter(corr_scores(:,4),corr_scores(:,5));
fit = lsline;
s.MarkerFaceColor = [0.5 0.5 0.5];
s.MarkerEdgeColor = [0 0 0];
s.MarkerFaceAlpha = 1;
s.MarkerEdgeAlpha = 1;
s.SizeData = 15;
fit.Color = [0 0 0];
fit.LineWidth = 2;
xlim([-3 5]); ylim([-2 2]);
set(gca,'xtick',[-2 0 2 4],'ytick',[-2 -1 0 1 2])
[roh(1,2),pval(1,2)] = corr(corr_scores(:,4),corr_scores(:,5),'tail','right');
text(1,-1.3,['coef = ',num2str(round(roh(1,2),3))],'color','k');
text(1,-1.5,['pval = ',num2str(round(pval(1,2),3))],'color','k');
xlabel('NC-C Naming (z-scores)'); ylabel('NC-C Semantic (z-scores)');
title('Reaction Times');
set(fig3,'position',[488.2 341.8 711.2 420],'color','w');

[Fisher(1,3),~,~,~] = corr_rtest(roh(1,1),roh(1,2),18,18);

%%
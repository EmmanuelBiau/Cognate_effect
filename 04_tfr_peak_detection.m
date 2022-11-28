%% Correct power spectrum and compute peak of power in theta, alpha and beta ranges;
clearvars;clc;
tfr_dir = 'XXX';
peak_out = 'XXX';

% find roi channels for naming and semantic tasks;
load stats_beta.mat; loal all_labels.mat;
chn_n = ismember(all_labels,stats_beta.L1N.roi);
chn_s = ismember(all_labels,stats_beta.L1S.roi);
clear stats_beta

subj_list = dir([tfr_dir,'*.mat']);
data_list = dir([tfr_dir,'*.mat']);

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

    % cycle through participants;
    for s = idx 

    name = subj_list(s).name;
    load([tfr_dir,name]);  

    % select only hit trials;
    cfg = [];
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.hit,1), freq_raw.trialinfo, 'UniformOutput', false)));
    freq_raw = ft_selectdata(cfg,freq_raw);

    % detect peaks of power;
    freq{1,1} = freq_raw;
    alltrials_spctrm = {};
    for ff = 1

        if strcmp(task,'L1N')==1
            chn = chn_n;
            t1 = 0.1;
            t2 = 0.22;            
            dist1 = abs(freq{1,1}.time - t1); dist2 = abs(freq{1,1}.time - t2);
            mindist1 = min(dist1); mindist2 = min(dist2);
            t1 = find(dist1 == mindist1); t2 = find(dist2 == mindist2);
            clear dist1 dist2 mindist1 mindist2
            
        elseif strcmp(task,'L1S')==1
            chn = chn_s;
            t1 = 0.28;
            t2 = 0.4;           
            dist1 = abs(freq{1,1}.time - t1); dist2 = abs(freq{1,1}.time - t2);
            mindist1 = min(dist1); mindist2 = min(dist2);
            t1 = find(dist1 == mindist1); t2 = find(dist2 == mindist2);
            clear dist1 dist2 mindist1 mindist2           
        end

            for trl = 1:size(freq{ff}.trialinfo,1)

                % select time-window of interest to compute the pow spectrum after stimulus onset;
                pow = squeeze(nanmean(nanmean(nanmean(freq{ff}.powspctrm(trl,chn,:,t1:t2),4),2),1));
                pow_temp = struct('time',1,...
                           'freq',freq{ff}.freq(freq{ff}.freq>=2&freq{ff}.freq<=40),...
                           'label',{{'dummy'}},...
                           'dimord','chan_freq_time',...
                           'powspctrm',pow(freq{ff}.freq>=2&freq{ff}.freq<=40)');

                % 1/f correction on rawpow; 
                pow_temp = uni_subtract1f(pow_temp);
                % detection of power peak in the spectrum ranges of interest (theta, alpha and beta);
                theta_peak.theta_peak(trl,1) = uni_getPeak(pow_temp.powspctrm',pow_temp.freq',[3 7],0);
                theta_peak.trialinfo = freq{ff}.trialinfo;

                alpha_peak.alpha_peak(trl,1) = uni_getPeak(pow_temp.powspctrm',pow_temp.freq',[8 12],0);
                alpha_peak.trialinfo = freq{ff}.trialinfo;            

                beta_peak.beta_peak(trl,1) = uni_getPeak(pow_temp.powspctrm',pow_temp.freq',[25 35],0);
                beta_peak.trialinfo = freq{ff}.trialinfo;

                % save corrected spectrum of all trials ;
                alltrials_spctrm{ff}.pow(trl,:) = pow_temp.powspctrm;
                alltrials_spctrm{ff}.onefcorr(trl,:) = pow_temp.onefcorr; 

            end

                alltrials_spctrm{ff}.freq = pow_temp.freq;
                alltrials_spctrm{ff}.label{1,1} = 'spect';
                alltrials_spctrm{ff}.dimord = 'rpt_freq';
                alltrials_spctrm{ff}.trialinfo = freq{ff}.trialinfo;
                alltrials_spctrm{ff}.time = 1;
                alltrials_spctrm{ff}.fractal = pow_temp.fractal; 

    end

    peak_info.theta = theta_peak;
    peak_info.alpha = alpha_peak;
    peak_info.beta = beta_peak;
    peak_info.correc_spectm = alltrials_spctrm;

    % fprintf('\nsave peaks of the subject...\n')
    save([peak_out,name(1:4),name(5:end-11),'_peak_info'],'peak_info','-v7.3');

    end

end

fprintf('\nPeak detection for all subjects done...\n')

% Find mean peaks L1N and L1S of each participant;
clear idx_subj
for ii = 1:18
    
    idx_subj(ii,1)= ii;
    sub = allsubj.id(cell2mat(allsubj.subj_num)==ii);
    idx_subj(ii,2) = sub{1,1}; idx_subj(ii,3) = sub{2,1};

end

allsubj_meanpeak_info = [];
for s = 1 : length(idx_subj)
    
    name1 = ['L1N_',num2str(idx_subj(s,2)),'_peak_info.mat']; load([peak_dir,name1]);
    peak_trials_theta = peak_info.theta.theta_peak;
    peak_trials_alpha = peak_info.alpha.alpha_peak;
    peak_trials_beta = peak_info.beta.beta_peak;
    clear peak_info
    name2 = ['L1S_',num2str(idx_subj(s,3)),'_peak_info.mat']; load([peak_dir,name2]);
    peak_trials_theta = [peak_trials_theta;peak_info.theta.theta_peak];
    peak_trials_alpha = [peak_trials_alpha;peak_info.alpha.alpha_peak];
    peak_trials_beta = [peak_trials_beta;peak_info.beta.beta_peak];
    clear peak_info
    allsubj_meanpeak_info(s,1) = s;
    allsubj_meanpeak_info(s,2) = mean(peak_trials_theta);
    allsubj_meanpeak_info(s,3) = mean(peak_trials_alpha);
    allsubj_meanpeak_info(s,4) = mean(peak_trials_beta);

end

% save([peak_dir,'allsubj_meanpeak_info'],'allsubj_meanpeak_info');

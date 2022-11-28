%% Analyse behavior results;
clearvars;clc; 
log_dir = 'XXX';
behav_out = 'XXX';

subj_list = dir([log_dir,'*.mat']);

allsubj = cell(length(subj_list),8);
allsubj = array2table(allsubj);
allsubj.Properties.VariableNames = {'subj_num' 'task' 'id' 'name' 'c_cr' 'c_rt' 'nc_cr' 'nc_rt'};

% cycle over tasks;
expe = {'L1N','L1S'}; task_num = 0;
for item = 1 : length(expe)
   
% select only pac from expe(item);
task = expe{item};
clear idx
for ii = 1 : length(subj_list)
    idx(1,ii) = contains(subj_list(ii).name,task);
end
idx = find(idx==1);

% cycle over participants;
count = 0;
for s = idx
  
count = count+1; 

% load data_clean and peak info of the participant;
name = subj_list(s).name;    
load([log_dir,name]); 

% compute cr and rt of the task;
c_cr = nanmean(logfile.cr(logfile.cognate==1));
c_rt = nanmean(logfile.rt(logfile.cognate==1 & logfile.cr==1 & logfile.rt < nanmean(logfile.rt)+2*nanstd(logfile.rt)...
         & logfile.rt > nanmean(logfile.rt)-2*nanstd(logfile.rt)));
     
nc_cr = nanmean(logfile.cr(logfile.cognate==2));
nc_rt = nanmean(logfile.rt(logfile.cognate==2 & logfile.cr==1 & logfile.rt < nanmean(logfile.rt)+2*nanstd(logfile.rt)...
         & logfile.rt > nanmean(logfile.rt)-2*nanstd(logfile.rt)));

allsubj.subj_num{count+task_num} = count;
allsubj.task{count+task_num} = task;
allsubj.id{count+task_num} = str2num(logfile.erp_code{1}(5:end));
allsubj.name{count+task_num} = logfile.name(1);
allsubj.c_cr{count+task_num} = c_cr;
allsubj.c_rt{count+task_num} = c_rt;
allsubj.nc_cr{count+task_num} = nc_cr;
allsubj.nc_rt{count+task_num} = nc_rt;

end
task_num = task_num+18;
end

% export behavior data;
save([behav_out,'allsubj_behav'],'allsubj');

% Calculate correlations accuracy difference (nc-c) between Picture-naming and Size-judgement tasks (to modify accordingly for the reaction times as well);
clear corr_scores 
corr_scores(:,1) = 1:18;
corr_scores(:,2) = (cell2mat(allsubj.nc_cr(strcmp(allsubj.task,'L1N')==1)) - cell2mat(allsubj.c_cr(strcmp(allsubj.task,'L1N')==1)))./cell2mat(allsubj.c_cr(strcmp(allsubj.task,'L1N')==1));
corr_scores(:,3) = (cell2mat(allsubj.nc_cr(strcmp(allsubj.task,'L1S')==1)) - cell2mat(allsubj.c_cr(strcmp(allsubj.task,'L1S')==1)))./cell2mat(allsubj.c_cr(strcmp(allsubj.task,'L1S')==1));
[roh,pval] = corrcoef(corr_scores(:,2),corr_scores(:,3));

close all; figure;
s = scatter(corr_scores(:,2),corr_scores(:,3));
fit = lsline;
s.MarkerFaceColor = [0.5 0.5 0.5];
s.MarkerEdgeColor = [0 0 0];
s.MarkerFaceAlpha = 1;
s.MarkerEdgeAlpha = 1;
s.SizeData = 15;
fit.Color = [0 0 0];
fit.LineWidth = 2;
xlim([-0.40 0.20]); ylim([-0.30 0.3]);
set(gca,'xtick',[-0.4 -0.2 0 0.2],'ytick',[-0.3 -0.15 0 0.15 0.3])
xlabel('Cognate-effect (L1N)')
ylabel('Cognate-effect (L1S)');
text(-0.38,0.2,['coef = ',num2str(round(roh(1,2),3))],'color','k');
text(-0.38,0.18,['pval = ',num2str(round(pval(1,2),3))],'color','k');

%%
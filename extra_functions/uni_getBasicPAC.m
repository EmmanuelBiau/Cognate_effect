function pac = uni_getBasicPAC(config,data)

% This function using the Modulation Index (Tort et al., 2010) to calculate
% phase-amplitude coupling between specified frequency bands. 

%% Prepare Data
% extract variables
if ~isfield(config,'latency'); latency = 'all';         else; latency = config.latency;     end
if ~isfield(config,'trials'); trials = 'all';           else; trials = config.trials;     end
if ~isfield(config,'phasefreq'); phasefreq = [];     else; phasefreq = config.phasefreq; end
if ~isfield(config,'powerfreq'); powerfreq = [];  else; powerfreq = config.powerfreq; end
if ~isfield(config,'nbins'); nbins = 18;                else; nbins = config.nbins; end
if ~isfield(config,'keeptrials'); keeptrials = 'no';    else; keeptrials = config.keeptrials; end

% get sizes
nchan = numel(data.label);
nphase = size(phasefreq,1);
npower = size(powerfreq,1);
nfreq = nphase + npower;

% select trials if requested
if strcmpi(trials,'all')
    ntrl = numel(data.trial);
else
    data = ft_selectdata(struct('trials',trials),data);
    ntrl = numel(data.trial);
end

% select time if required
if strcmpi(latency,'all')
    ntime = numel(data.time{1});
else
    % FIX ME: figure out how fieldtrip selects latency because it's not
    % the commented out method:    
    %ntime = sum(round(data.time{1},4)>=latency(1) & round(data.time{1},4)<=latency(end));
    
    cfg = [];
    cfg.latency = latency;
    tmp = ft_selectdata(cfg,data);
    ntime = numel(tmp.time{1});
end

%% Filter Data
% define filter bands
combi_band = cat(1,phasefreq,powerfreq);

% predefine signal
signal = zeros(ntrl,nchan,ntime,size(combi_band,1),'single');

% cycle through each band frequency
fprintf('\nfiltering data...\n') 
parfor freq = 1 : nfreq
    
    % extract frequency band
    tmp_freq = combi_band(freq,:);
    
    % filter
    cfg             = [];
    cfg.hpfilter    = 'yes';
    cfg.hpfreq      = tmp_freq(1);
    cfg.hpfilttype  = 'but';
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = tmp_freq(2);
    cfg.lpfilttype  = 'but';
    filtered        = ft_preprocessing(cfg,data);

    % select latency
    cfg             = [];
    cfg.latency     = latency;
    filtered        = ft_selectdata(cfg,filtered);

    % timelock gamma data
    cfg             = [];
    cfg.keeptrials  = 'yes';
    filtered        = ft_timelockanalysis(cfg,filtered);
    signal(:,:,:,freq) = single(filtered.trial);    
end

%% Get Power/Phase
% predefine output measures
phapow = zeros(ntrl,nchan,nphase);
powpow = zeros(ntrl,nchan,npower);

% get theta phase/gamma amp.
fprintf('\nextracting phase and amplitude metrics...\n');
for trl = 1 : ntrl
    for chan = 1 : nchan
        
        % get output power
        phapow(trl,chan,:) = mean(abs(hilbert(squeeze(signal(trl,chan,:,1:nphase)))),1);
        powpow(trl,chan,:) = mean(abs(hilbert(squeeze(signal(trl,chan,:,nphase+1:end)))),1);       
        
        % get phase/power for analysis
        signal(trl,chan,:,1:nphase) = angle(hilbert(squeeze(signal(trl,chan,:,1:nphase))));
        signal(trl,chan,:,nphase+1:end) = abs(hilbert(squeeze(signal(trl,chan,:,nphase+1:end))));
    end
end

%% Calculate MI
% extract trials of interest
fprintf('binning signal amplitude...\n');
[mi,bin_norm,~,edges] = int_calculateMI(signal,nbins,nchan,nphase,npower,ntrl,keeptrials);

%% Consolidate
% add data into pac structure
fprintf('packaging data...\n')
pac             = [];
pac.label       = data.label;
pac.freq        = mean(powerfreq,2);
pac.time        = mean(phasefreq,2);
pac.cfg         = [];
pac.dimord      = 'chan_freq_time';
pac.powspctrm   = permute(mi,[2 4 3 1]);
pac.phase       = bin_norm;
pac.powpow      = powpow;
pac.phapow      = phapow;
pac.bin_val     = edges;
pac.bin_val     = pac.bin_val(2:end) - diff(pac.bin_val)/2;
if strcmpi(keeptrials,'yes'); pac.trialinfo = data.trialinfo; end

end

%% Subfunctions
function [bin_norm,edges] = int_binPower(signal,bin_edge)

% get theta bins
binned_amp = zeros(bin_edge,1);

% get bin indices
[~,edges,bin_val] = histcounts(signal(:,1),bin_edge,'BinLimits',[-pi,pi]);

% cycle through each bin
t = unique(bin_val);
t(t==0) = [];
for i = 1:numel(t)
    
    % get mean amplitude
    binned_amp(t(i),1) = sum(signal(bin_val==t(i),2)) ./ sum(bin_val==t(i));    
end

% normalise across bins
bin_norm = binned_amp ./ sum(binned_amp);
end

function [mi,bin_norm,bin_perm,edges] = int_calculateMI(signal,nbins,nchan,nphase,npower,ntrl,keeptrials)

if strcmpi(keeptrials,'no')

    % predefine MI matrix       
    bin_norm    = zeros(nbins,nchan,nphase,npower);
    bin_perm    = zeros(100,nbins,nchan,nphase,npower);

    % cycle through each trial
    tic
    for chan = 1 : nchan

        % restrict signal to current channel
        sigtmp = permute(signal(:,chan,:,:),[1 3 4 2]);

        % collapse over trials
        sigtmp = permute(sigtmp,[3 1 2]);
        sigtmp = sigtmp(:,:);

        % cycle through each specified phase and power
        for i = 1 : nphase
            for j = 1 : npower

                % extract signal of interest
                soi = sigtmp([i j+nphase],:)';

                % bin power
                [bin_norm(:,chan,i,j),edges] = int_binPower(soi,nbins);
            end
        end 
        
        if chan == 10
            te = toc;
            tr = (te/10)*(nchan-10);
            fprintf('estimated time remaining: %d minutes...\n',round(tr/60))
        end
    end

    % tidy
    clear samp bin_val binned_amp bin_quant sig_tmp

    % get modulation index for observed data
    fprintf('calculating MI...\n');
    p = bin_norm;
    q = zeros(size(p)) + repmat(mean(p,1),[size(p,1) 1 1 1]);
    mi = (sum(p .* log10(p ./ q),1)) ./ log10(size(p,1));
    
else
    
    % predefine MI matrix       
    bin_norm    = zeros(nbins,ntrl,nchan,nphase,npower);
    bin_perm    = zeros(100,ntrl,nbins,nchan,nphase,npower);

    % cycle through each trial
    tic
    for trl = 1 : ntrl
        for chan = 1 : nchan

            % restrict signal to current channel
            sigtmp = squeeze(signal(trl,chan,:,:));

            % cycle through each specified phase and power
            for i = 1 : nphase
                for j = 1 : npower

                    % extract signal of interest
                    soi = sigtmp(:,[i j+nphase]);

                    % bin power
                    [bin_norm(:,trl,chan,i,j),edges] = int_binPower(soi,nbins);
                end
            end 
        end
        if trl == 10
            te = toc;
            tr = (te/10)*(ntrl-10);
            fprintf('estimated time remaining: %d minutes...\n',round(tr/60))
        end
    end

    % tidy
    clear samp bin_val binned_amp bin_quant sig_tmp

    % get modulation index for observed data
    fprintf('calculating MI...\n');
    p = bin_norm;
    q = zeros(size(p)) + repmat(mean(p,1),[size(p,1) 1 1 1 1]);
    mi = permute(((sum(p .* log10(p ./ q),1)) ./ log10(size(p,1))),[2 3 4 5 1]);
        
end
end

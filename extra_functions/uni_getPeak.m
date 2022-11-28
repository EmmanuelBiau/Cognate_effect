function peak = uni_getPeak(pow,freq,frange,detrend)

% restrict to frange
pow     = pow(freq>=frange(1)&freq<=frange(2));
freq    = freq(freq>=frange(1)&freq<=frange(2));

% find theta peak
[val,idx] = findpeaks(pow);
rep       = 1;

% if detrend first
if exist('detrend','var'); rep = detrend; idx=[]; end

% detrend if no peak found
while isempty(idx)

    % detrend
    pft         = polyfit(freq,pow,rep);
    pval        = polyval(pft,freq);
    diff        = pow - pval;

    % find peak
    [~,idx,~,val] = findpeaks(diff);
    rep = rep + 1;

    % break if rep > 3
    if rep > 3; break; end
end

% get max peak
[~,max_idx] = max(val);
peak = freq(idx(max_idx));
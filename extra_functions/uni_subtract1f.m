function freq = uni_subtract1f(freq)

% This function estimates and subtracts the fractal property of a power 
% spectrum. The [freq] input must match a Fieldtrip formatted
% time-frequency data structure.

% extract parameters
f   = freq.freq;
pow = freq.powspctrm;

% exclude bandstop filtered frequencies from fit
fidx = f>45&f<55 | f>95;
f(fidx) = [];
pow(fidx) = [];
freq.freq(fidx) = [];

% fit 1/f
logf    = squeeze(log10(f));
logpow  = squeeze(log10(pow));
beta    = zeros(size(logpow,1),1);
onefcorr = zeros(size(logpow,1),length(log10(f)));
% cycle through each trial
for trl = 1 : size(logpow,1)
   
    % get temporary power and frequency
    tmp_f   = logf;
    tmp_pow = logpow(trl,:);
    
    % iterate
    while true
    
        % get fit
        b = [ones(size(tmp_f))' tmp_f'] \ tmp_pow';
        linft = tmp_f*b(2) + b(1);

        % subtract fit
        firstpass = tmp_pow - linft;

        % get error (the mean power of all frequencies less than zero)
        erange = mean(firstpass(firstpass<0));
        
        % find all postive values greater than the error
        eidx   = firstpass>abs(erange);
        
        % if no error
        if sum(eidx) == 0
        
            % recompute fit using all data
            linft = logf*b(2) + b(1);
            
            % subtract fit
            logpow(trl,:) = logpow(trl,:) - linft;
                        
            % save output 
            beta(trl,1) = b(2);
            onefcorr(trl,:) = linft;
            break
            
        % else if more than half of the frequencies have been removed
        elseif numel(tmp_f) <= numel(logf)/2
            
            % recompute fit using all data
            linft = logf*b(2) + b(1);
            
            % subtract fit
            logpow(trl,:) = logpow(trl,:) - linft;
                        
            % save output 
            beta(trl,1) = b(2);
            onefcorr(trl,:) = linft;
            break
                        
         % otherwise, remove freqencies that exceed the error threshold 
        else
            tmp_pow = tmp_pow(eidx==0);
            tmp_f   = tmp_f(eidx==0);
        end
    end
end

% convert pow
freq.powspctrm = 10.^logpow;
freq.fractal = repmat(beta,[1 size(freq.powspctrm,2)]);
freq.onefcorr = onefcorr;

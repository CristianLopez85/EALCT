function [tfc, tfrtic, tcrtic] = sqSTCTm(x, resol, tDS, g)
% x: input signal
% resol: resolution for frequency and chirp rate
% tDS: hopping in time
% h: window function, a column vector

[xrow,xcol] = size(x) ;
t = 1:length(x) ;
tLen = length(t(1:tDS:length(x))) ;

	% for tfr
M = length(-0.5+resol:resol:0.5) ;
crate = ((1:M-1)-ceil(M/2))/M^2; % discretization of chirp rate

	% for tfrsq
cLen = length(crate); % number of bins of chirp rate
%====================================================================
	%% check input signals
if (xcol~=1)
    error('X must have only one column');
elseif (tDS < 1) || (rem(tDS,1)) 
    error('tDS must be an integer value >= 1');
end

[hrow,hcol] = size(g); Lh = (hrow-1)/2; 
if (hcol~=1)||(rem(hrow,2)==0)
    error('H must be a smoothing window with odd length');
end
ht = -Lh:Lh ;

%====================================================================
	%% run STFT and reassignment rule
tfc = zeros(M-1, M/2, tLen); 	% chirplet transform

% lowFreq: minimal frequency to observe, smallest as 0
% highFreq: maximal frequency to observe, largest as 0.5

tfrtic = linspace(0, 0.5, M/2)' ; % positive frequency 
tcrtic = crate;

% fprintf(['Chirp-rate total: ',num2str(cLen), '; now:     ']) ;            Commented for rebbutal

for cidx = 1:M-1
%     fprintf('\b\b\b\b') ;	tmp = sprintf('%4d',cidx) ; fprintf(tmp) ;      Commented for rebbutal
    chirp = crate(cidx);
    for tidx = 1:tLen
        
        % ti is the current time
        ti = t((tidx-1)*tDS+1);
        % tau is the relevant index associated with ti
        tau = -min([round(M/2)-1,Lh,ti-1]):min([round(M/2)-1,Lh,xrow-ti]);
        % indices is the absolute index in the "evaluation window"
        indices= rem(M+tau,M)+1;
        tf0 = zeros(M, 1) ;
        tf0(indices) = x(ti+tau).*conj(g(Lh+1+tau)).*exp(-pi*1i*chirp.*(ht(Lh+1+tau)').^2);
        % select nonnegative frequencies
        tf0 = fft(tf0) ; % tf0 = tf0(1:M/2) ;
        tfc(cidx, :, tidx) = tf0(1:M/2) ;
        
    end
end
% fprintf('\n') ;                                                          Commented for rebbutal
end
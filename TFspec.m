function [ASpec fbin] = TFspec(IFmulti,IAmulti,band)
% Construct an adaptive time-frequency spectrum (ATFS) using the estimated IAs and IFs
% band = [min max] denotes the frequency range (from min to max) shown in the ATFS 
% IFmulti: estimated instantaneous frequency (IF) time series for all the signal modes; each row of IFmulti corresponds to the IF of each mode
% IAmulti: estimated instantaneous amplitude (IA) time series for all the signal modes; each row of IAmulti corresponds to the IA of each mode
frnum = 1024; %the number of the frequency bins 
fbin = linspace(band(1),band(2),frnum);
num = size(IFmulti,1); % the number of the signal modes
N = size(IFmulti,2); %length of the signal
ASpec = zeros(frnum,N);
delta = floor(frnum*0.1e-2);%
for kk = 1:num
    temp = zeros(frnum,N);
    for ii = 1:N
        [~,index] = min(abs(fbin - IFmulti(kk,ii)));
        lindex = max(index-delta,1);
        rindex = min(index+delta,frnum);
        temp(lindex:rindex,ii) = IAmulti(kk,ii);
    end
    ASpec = ASpec + temp;
end
        
    



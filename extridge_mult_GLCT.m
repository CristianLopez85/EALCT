function [fidexmult, tfdv] = extridge_mult_GLCT(Sig, SampFreq, num, delta, orderIF,bw,alpha,N,hlength)
% Extract ridges for multi-component signals.
% In each iteration,the signal component associated with the extrated ridge is
% reconstructed by the ICCD and then removed from the original signal so
% that the ridge curves of other signal components with smaller energies
% can be extracted in the subsequent iterations.
%%%%%%%%%%%%%%%%%%%%%%%    input      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sig��measured signal,a row vector
% SampFreq: sampling frequency
% num: the number of the signal components
% delta��maximum allowable frequency variation between two consecutive points
% orderIF: the order of the Fourier model used for smoothing the extracted ridge curves
% bw��the bandwidth of the ICCD (unit��Hz); herein the ICCD can be regarded as a time-frequency filtering technique
% Nfrebin,window are two parameters for implementing the STFT
% alpha��Tikhonov regularization parameter for ICCD.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  output   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fidexmult: obtained ridge index for multi-component signals
% tfdv: the corresponding ridge magnitude 

if (isreal(Sig))
Sig = hilbert(Sig);
end

orderamp = round(bw*length(Sig)/SampFreq);%Fourier order for characterizing signal amplitude
fidexmult = zeros(num,length(Sig));
tfdv = zeros(num,length(Sig));
for i = 1:num
% [Spec,f] = STFT(Sig(:),SampFreq,Nfrebin,window); %STFT

%--------------------------------------------------------------------------
[Spec,~,f] = GLCT((real(Sig(:))),N,SampFreq,hlength);
% Spec = Spec';
% [nLevel, ~] = size(Spec);
% f = linspace(-SampFreq/2,SampFreq/2,nLevel);
%--------------------------------------------------------------------------

c = findridges(Spec,delta);
[extr_Sig,~,IFfit] = ICCD(Sig,SampFreq,f(c),orderIF,orderamp,alpha);
findex = zeros(1,length(Sig));
for j = 1:length(Sig)
   [~,findex(j)] = min(abs(f - IFfit(1,j)));
   tfdv(i,j) = abs(Spec(findex(j),j));
end
fidexmult(i,:) = findex;
Sig = Sig - extr_Sig; % remove the extracted signal component so that other ridge curves with smaller energies can be extracted in the subsequent iterations
end
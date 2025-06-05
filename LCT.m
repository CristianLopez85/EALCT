function [tfr] = LCT(x,c,fs,h)
%      Linear Chirplet Transform
%	x     : Signal.
%	c     : Chirp Rate.
%	fs    : Sample Frequency.
%	h     : Window Function.

%	tfr   : Time-Frequency Representation.

%  This program is free software; you can redistribute it and/or modify
%  it according to your requirement.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
%   Written by Gang Yu in Shandong University at 2015.4.28.

[xrow,~] = size(x);

if (nargin < 3)
    error('At least 3 parameter is required');
end

N=xrow;
t=1:xrow;

[~,tcol] = size(t);

[hrow,~]=size(h); Lh=(hrow-1)/2;

tt=(1:N)/fs;
tfr= zeros (N,tcol) ;

for icol=1:tcol
    ti= t(icol); 
    tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
    indices= rem(N+tau,N)+1;

    rSig = x(ti+tau,1);    %rSig = Hilbert(real(rSig));

    tfr(indices,icol)=rSig.*conj(h(Lh+1+tau)).*exp(-1i * 2 * pi * (c/2) * (tt(ti+tau)-tt(icol)).^2)';

end
tfr=fft(tfr);
tfr=tfr(1:round(end/2),:);
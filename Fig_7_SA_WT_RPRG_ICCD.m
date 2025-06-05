clear;clc;close all
SampFreq = 400; tend = 2;
t = 1/SampFreq : 1/SampFreq : tend;

Sig1 = exp(2*pi*1i*(-20*t.^2 + 90*t));               % IF1 = -40*t + 90;
Sig2 = exp(2*pi*1i*(20*t.^2 + 10*t));                % IF2 =  40*t + 10;

% Sig = Sig1 + Sig2;
%
Ratio = 0;
% window = 128;
% Nfrebin = 1024;
%
bw1 = 0.01*SampFreq; % SampFreq/60;%
orderIF1 = 50; % use the Fourier model to smooth the detected ridge curves£»orderIF1 could be larger than the following orderIF
num = 2; % the number of the components
delta = 0.03*SampFreq; % 20; % ;                                            % attempt by CL 0.03*SampFreq
alpha1 = 5;
%
alpha2 = 0.5;
bw = SampFreq/40; %  % bandwidth of the ICCD
orderIF = 5;

j = 0; pin = 10;        % pin = 9: fails, so manually start at 10.
for SNR = 15%-4:1:15
    j = j + 1;
    for k = 1 : 10
        pin = pin + 1;
        Sig = Sig1 + Sig2;
        seed = pin; % Use the loop variable as the seed
        rng(seed); % Set the seed for random number generation
        Sig = awgn(Sig,SNR);

        %%%%%%%%%%%STFT%%%%%%%%%%%%%%%%%%%%%
        N = 7; hlength = 100;
        [tfr,~,f] = GLCT((real(Sig(:))),N,SampFreq,hlength);
        [fr,mods] = WaveletSignal(t,real(Sig(:)),0,100,100,2,'e','e');
        mods = mods/max(mods,[],'All');

        subplot(3,1,2)
        MoDAL.WTSpectraPlot(time,freq,mods,options)
        %         [nLevel, ~] = size(tfr);
        %         f = linspace(-SampFreq/2,SampFreq/2,nLevel);
        %         [Spec,f] = STFT(Sig(:),SampFreq,Nfrebin,window);

        % ridge extraction and fitting
        [fidexmult, tfdv] = extridge_mult_GLCT(Sig, SampFreq, num, delta, orderIF1,bw1,alpha1,N,hlength);

        % ridge path regrouping
        thrf = 0.03*SampFreq; % length(f)/30; %                            %  attempt by by CL  0.03*SampFreq

        [findex,interset] = RPRG(fidexmult,thrf);

        % signal decomposition
        orderamp = round(bw*length(Sig)/SampFreq);% Fourier order for characterizing signal amplitudes
        [extr_Sig,ampmatrix,IFfit] = ICCD(Sig,SampFreq,f(findex),orderIF,orderamp,alpha2); % ICCD performs the joint-optimization scheme using all the obtained IFs

        % by CL
        Sig1_s = Sig1(100: end-100); Sig2_s = Sig2(100: end-100);
        mode_s = extr_Sig(1,100: end-100); mode2_s = extr_Sig(2,100: end-100);

        e1 = norm(real(Sig1_s-mode_s))/norm(real(mode_s));
        e2 = norm(real(Sig2_s-mode2_s))/norm(real(mode2_s));

        e1_r = norm(real(Sig2_s-mode_s))/norm(real(mode_s));               % Because, may Sig1_s can be displayed as Sig2_s
        e2_r = norm(real(Sig1_s-mode2_s))/norm(real(mode2_s));

        e_k_1 = mean([e1 e2]);
        e_k_2 = mean([e1_r e2_r]);

        e_k(k) = min([e_k_1 e_k_2]);
    end
    e(j) = mean(e_k);
    j
end
e=e';

function [fr,mods] = WaveletSignal(time,signal,minFreq,maxFreq,numFreq,motherWaveletFreq,mirrori,mirrorf)
chp1 = 0.2;
chp2 = 0.2;

[x_mirror,L_chp1,L_chp2,NoMirrorIni,NoMirrorEnd] = MirrorSignal(time,signal,mirrori,mirrorf,chp1,chp2);

[fr,mods] = MoDAL.WaveletTransform(time,x_mirror,minFreq,maxFreq, ...
    'motherWaveletFreq',motherWaveletFreq,'numFreq',numFreq);

if NoMirrorIni == 0 && NoMirrorEnd == 0
    mods = mods(L_chp1:end-L_chp2+1,:);
elseif NoMirrorIni == 0 && NoMirrorEnd == 1
    mods = mods(L_chp1:end);
elseif NoMirrorIni == 1 && NoMirrorEnd == 0
    mods = mods(1:end-L_chp2+1);
end
end

function [x_mirror,L_chp1,L_chp2,NoMirrorIni,NoMirrorEnd] = ...
    MirrorSignal(time,signal,mirrori,mirrorf,chp1,chp2)
arguments
    time (:,1) double
    signal (:,1) double
    mirrori string = 'e';
    mirrorf string  = 'e';
    chp1 double = 0.2;
    chp2 double = 0.2;
end
L = length(time); L_chp1 = chop(L*chp1,2); L_chp2 = chop(L*chp2,2);
[x_o_ini,~] = MirrorImgSigOIni(time(1:L_chp1),signal(1:L_chp1));
[x_e_ini,~] = MirrorImgSigEIni(time(1:L_chp1),signal(1:L_chp1));
[x_o_fin,~] = MirrorImgSigOFin(time(end-L_chp2+1:end),signal(end-L_chp2+1:end));
[x_e_fin,~] = MirrorImgSigEFin(time(end-L_chp2+1:end),signal(end-L_chp2+1:end));

NoMirrorIni = 0;
NoMirrorEnd = 0;
if strcmp(mirrori,'e')
    if strcmp(mirrorf,'e')
        x_mirror = [x_e_ini; signal; x_e_fin];
    elseif strcmp(mirrorf,'o')
        x_mirror = [x_e_ini; signal; x_o_fin];
    else
        x_mirror = [x_e_ini; signal];
        NoMirrorEnd = 1;
    end
elseif strcmp(mirrori,'o')
    if strcmp(mirrorf,'e')
        x_mirror = [x_o_ini; signal; x_e_fin];
    elseif strcmp(mirrorf,'o')
        x_mirror = [x_o_ini; signal; x_o_fin];
    else
        x_mirror = [x_o_ini; signal];
        NoMirrorEnd = 1;
    end
else
    if strcmp(mirrorf,'e')
        x_mirror = [signal; x_e_fin];
    elseif strcmp(mirrorf,'o')
        x_mirror = [signal; x_o_fin];
    else
        x_mirror = signal;
        NoMirrorEnd = 1;
    end
    NoMirrorIni = 1;
end
end

function X = chop(Xin,n,unit)
%CHOP   CHOP(X,n) rounds elements of X to n significant figures.
%       CHOP(X,n,unit) rounds the elements of X to n significant
%   figures whose digits (mantissa) are exactly divisible
%   by unit.
%
%       e.g. chop(3.141592,5)   returns 3.141600000..
%       e.g. chop(3.141592,3,5) returns 3.150000000..
%            chop(3.141592,3,3) returns 3.150000000..
%            chop(3.141592,3,2) returns 3.140000000..
%

%   Copyright 1986-2002 The MathWorks, Inc.

% Set last sig. fig. rounding to 1 if only two input arguments.
if nargin<3
    unit=1;
end

% Cater for -ve numbers  and numbers = 0.
X = abs(Xin) +(Xin==0);
[nx,mx] = size(X);
exponent = unit.*((10*ones(nx,mx)).^(floor(log10(X))-n+1));
X = round(X./exponent).*exponent;

% Put back sign and zeros
X = sign(Xin).*X.*(Xin~=0);

end

function [x_mirror, t_mirror] = MirrorImgSigOIni(t, x)
% Generate the mirror image signal by an odd symmetry about t=t0
L = length(t); xtmp = x - x(1);
t_mirror = [-(t(L:-1:2)-t(1))+t(1)];
x_mirror = [-xtmp(L:-1:2)+x(1)];
end

function [x_mirror, t_mirror] = MirrorImgSigEIni(t, x)
% Generate the mirror image signal by an even symmetry about t=t0
L = length(t); xtmp = x - x(1);
t_mirror = [-(t(L:-1:2)-t(1))+t(1)];
x_mirror = [xtmp(L:-1:2)+x(1)];
end

function [x_mirror, t_mirror] = MirrorImgSigOFin(t, x)
% Generate the mirror image signal by an odd symmetry about t=tf
dt = t(2)-t(1); L = length(t); xtmp = x - x(end);
t_mirror = [t(end)+dt*[1:L-1]'];
x_mirror = [-xtmp(end-1:-1:end-L+1)+x(end)];
end

function [x_mirror, t_mirror] = MirrorImgSigEFin(t, x)
% Generate the mirror image signal by an even symmetry about t=tf
dt = t(2)-t(1); L = length(t); xtmp = x - x(end);
t_mirror = [t(end)+dt*[1:L-1]'];
x_mirror = [xtmp(end-1:-1:end-L+1)+x(end)];
end
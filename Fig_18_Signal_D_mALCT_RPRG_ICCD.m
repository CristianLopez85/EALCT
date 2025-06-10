clear;clc;close all

SampFreq = 400; tend = 2;
t = 1/SampFreq : 1/SampFreq : tend; Sf = 100;% time = time(1:400);
% 
amp1 = exp(0.1*t);  amp2 = exp(0.2*t);  amp3 = exp(0.3*t);
% 
Sig1 = amp1.*exp(2*1i*pi*(25*t.^2));                        IF1 = 50*t;
Sig2 = amp2.*exp(2*1i*pi*(25*t));                           IF2 = 25*ones(length(t),1);
Sig3 = amp3.*exp(2*1i*pi*(-22*t.^2 + 95*t + sin(5*t)));     IF3 = -44*t + 95 + 5*cos(5*t);

% Sig = Sig1 + Sig2;
%
Ratio = 0;
% window = 128;
% Nfrebin = 1024;
%
bw1 = SampFreq/60;
orderIF1 = 50; % use the Fourier model to smooth the detected ridge curves£»orderIF1 could be larger than the following orderIF
num = 2; % the number of the components
delta = 0.03*SampFreq; %20; % ;                                                           % attempt by CL 0.03*SampFreq
alpha1 = 5;
%
alpha2 = 0.5;
bw = SampFreq/40; %  % bandwidth of the ICCD       
orderIF = 5;

Nh = 151; Nc = 71;

j = 0; pin = 0; 
for SNR = -4:1:15
    j = j + 1;
    for k = 1 : 10
        pin = pin + 1;
        Sig = Sig1 + Sig2;
        seed = pin; % Use the loop variable as the seed
        rng(seed); % Set the seed for random number generation
         Sig = awgn(Sig,SNR);

        %%%%%%%%%%%STFT%%%%%%%%%%%%%%%%%%%%%
        %         N = 7; hlength = 100;
        %         tfr = GLCT((real(Sig(:))),N,SampFreq,hlength);

        Nx=length(Sig); Nfext=Nx;

        %Atau, Af and Ac are arguments of LCT time, frequency and chirp rate, respectively.
        Afext=(0:(Nfext-1))*SampFreq/Nfext;
        Nf=round(Nfext/2);
        f = Afext(1:Nf);
%         [Spec,~,f] = zALCT(Sig,SampFreq,Nh,Nc,num);
        %         [nLevel, ~] = size(Spec);
%         f = linspace(-SampFreq/2,SampFreq/2,nLevel);
%         [Spec,f] = STFT(Sig(:),SampFreq,Nfrebin,window);

        % ridge extraction and fitting
        [fidexmult, tfdv] = extridge_mult_ALCT(Sig, SampFreq, num, delta, orderIF1,bw1,alpha1,Nh,Nc);

        % ridge path regrouping
        thrf = 0.03*SampFreq; % length(f)/30; %                                              %  attempt by by CL  0.03*SampFreq 
        
        [findex,interset] = RPRG(fidexmult,thrf);
%         figure;subplot(211);plot(fidexmult(1,:));subplot(212);plot(findex(2,:))
        % signal decomposition
        orderamp = round(bw*length(Sig)/SampFreq);% Fourier order for characterizing signal amplitudes

        findex (findex < 1) = 1;                                           % Added for rebuttal

        [extr_Sig,ampmatrix,IFfit] = ICCD(Sig,SampFreq,f(findex),orderIF,orderamp,alpha2); % ICCD performs the joint-optimization scheme using all the obtained IFs

        % by CL
        chops = 100;
        Sig1_s = Sig1(chops: end-chops); Sig2_s = Sig2(chops: end-chops); Sig3_s = Sig3(chops: end-chops);
        mode_s = extr_Sig(1,chops: end-chops); mode2_s = extr_Sig(2,chops: end-chops); mode3_s = extr_Sig(2,chops: end-chops);

        e1 = norm(real(Sig1_s-mode_s))/norm(real(mode_s));
        e2 = norm(real(Sig2_s-mode2_s))/norm(real(mode2_s));
        e3 = norm(real(Sig3_s-mode3_s))/norm(real(mode3_s));

        e1_r = norm(real(Sig2_s-mode_s))/norm(real(mode_s));               % Because, may Sig1_s can be displayed as Sig2_s
        e2_r = norm(real(Sig3_s-mode2_s))/norm(real(mode2_s));
        e3_r = norm(real(Sig1_s-mode3_s))/norm(real(mode3_s));

        e1_r1 = norm(real(Sig3_s-mode_s))/norm(real(mode_s));               % Because, may Sig1_s can be displayed as Sig2_s
        e2_r1 = norm(real(Sig1_s-mode2_s))/norm(real(mode2_s));
        e3_r1 = norm(real(Sig2_s-mode3_s))/norm(real(mode3_s));
        
        e_M1 = [e1 e1_r e1_r1];
        e_M2 = [e2 e2_r e2_r1];
        e_M3 = [e3 e3_r e3_r1];

        e1_s = min(e_M1); e2_s = min(e_M2); e3_s = min(e_M3);

        e_k(k) = mean([e1_s e2_s e3_s]);
    end
    e(j) = mean(e_k);
    j
end
e=e';

clc; clear; close all
warning('off')
SampFreq = 100;
t = 0:1/SampFreq:3; 

Sig1 = exp(2*pi*1i*(30*t + 3*sin(3*t)));                                   %  IF2 = 30 + 9*cos(3*t);
Sig2 = 0.9*exp(2*pi*1i*(30*t - 3*sin(3*t) + 3/(2*pi)));                    %  IF1 = 30 - 9*cos(3*t);

% Sig = Sig1 + Sig2;
% rng('default')
%
bw1 = SampFreq/60;%
orderIF1 = 50; % use the Fourier model to smooth the detected ridge curves£»orderIF1 could be larger than the following orderIF
num = 2; % the number of the components
delta = 0.03*SampFreq; % ;                                                              % attempt by CL 0.03*SampFreq
alpha1 = 5;
%
alpha2 = 0.5;
bw = SampFreq/40; %  % bandwidth of the ICCD                         
orderIF = 5;

Nh = 71; Nc = 31;  % 151,71

j = 0; pin = 10; 
for SNR = -3:1:15
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
        thrf = 0.03*SampFreq; %length(f)/30; %                                         %  attempt by by CL  0.03*SampFreq 
        
        [findex,interset] = RPRG(fidexmult,thrf);
%         figure;subplot(211);plot(fidexmult(1,:));subplot(212);plot(findex(2,:))
        % signal decomposition
        orderamp = round(bw*length(Sig)/SampFreq);% Fourier order for characterizing signal amplitudes

        findex (findex < 1) = 1;                                           % Added for rebuttal

        [extr_Sig,ampmatrix,IFfit] = ICCD(Sig,SampFreq,f(findex),orderIF,orderamp,alpha2); % ICCD performs the joint-optimization scheme using all the obtained IFs

        % by CL
        Sig1_s = Sig1(25: end-25); Sig2_s = Sig2(25: end-25);
        mode_s = extr_Sig(1,25: end-25); mode2_s = extr_Sig(2,25: end-25);

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

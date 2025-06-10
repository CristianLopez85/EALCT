clear;clc;close all
SampFreq = 400; tend = 2;
t = 1/SampFreq : 1/SampFreq : tend;

Sig1 = exp(2*pi*1i*(-20*t.^2 + 90*t));               % IF1 = -40*t + 90;
Sig2 = exp(2*pi*1i*(20*t.^2 + 10*t));                % IF2 =  40*t + 10;

initstate(1) ;
am1 = smooth(cumsum(randn(length(t),1)) ./ SampFreq, 200, 'loess') ;
am1 = 2 + am1 ./ max(abs(am1)) ; am1 = am1';
am2 = smooth(cumsum(randn(length(t),1)) ./ SampFreq, 200, 'loess') ;
am2 = 2 + am2 ./ max(abs(am2)) ; am2 = am2';

rng('default')

Sig1 = am1.*Sig1;
Sig2 = am2.*Sig2;

% Sig = Sig1 + Sig2;

bw1 = SampFreq/60;%
orderIF1 = 50; % use the Fourier model to smooth the detected ridge curves£»orderIF1 could be larger than the following orderIF
num = 2; % the number of the components
delta = 0.03*SampFreq; % 20; % ;                                            % attempt by CL 0.03*SampFreq
alpha1 = 5;
%
alpha2 = 0.5;
bw = SampFreq/40; % bandwidth of the ICCD                            % attempt by CL 0.01*SampFreq
orderIF = 5;

[Siglength,~] = size(Sig1');
f = 0:SampFreq/Siglength:(SampFreq/2)-(SampFreq/2)/Siglength;  % From line 46

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
        N = 7; hlength = 100;
%         [tfr,~,f] = GLCT((real(Sig(:))),N,SampFreq,hlength);
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

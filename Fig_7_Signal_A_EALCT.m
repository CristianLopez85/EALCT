clc; clear; close all
warning('off')
tic
Hz = 400;
t = 0:1/Hz:2; t = t';

% figure; set(gcf,'Position',[20 100 800 450]);

Sig1 = exp(2*pi*1i*(-20*t.^2 + 90*t));                   IF1 = -40*t + 90;
Sig2 = 0.6*exp(2*pi*1i*(20*t.^2 + 10*t));                IF2 =  40*t + 10;

Nh = 151; Nc = 71;

%%
j = 0; pin = 0;
for SNR = -4%:1:15
    j = j + 1;
    for k = 1 : 10
        pin = pin + 1;
        Sig = Sig1 + Sig2;
        seed = pin; % Use the loop variable as the seed
        rng(seed); % Set the seed for random number generation
        Sig = awgn(Sig,SNR);

        % ridge extraction and fitting
        n = 2;
        Num = max(floor(n));    % the number of the components
        deltaf = 0.03*Hz; 
        beta = 1e-3;
        Sig0 = real(Sig)';
        if (isreal(Sig0))
            Sig0 = hilbert(Sig0);
        end

        fidexmult = zeros(Num,length(Sig0)); clearvars Spec Af c bw IFfit extr_Sig findex fidexmult
        num = Num;
        for i = 1 : Num

            [Spec,Atau,Af] = zALCT(Sig0,Hz,Nh,Nc,num);
            Spec = Spec';
            c = findridges_m(Spec,deltaf,0.5,num,7);

            c (c > round(length(t)/2)) = round(length(t)/2);                              % Added for rebuttal
            c (c < 1) = 1;                                                 % Added for rebuttal

            bw = 0.01*Hz; 
            [IFfit,extr_Sig] = Dechirp_filter(Sig0,Hz,bw,Af(c),beta);
            fidexmult(i,:) = c;
            Sig0 = Sig0 - extr_Sig;
            num = num - 1;
        end

        % ridge path regrouping
        Df = 0.03*Hz; 
        [findex,interset] = RPRG(fidexmult,Df);

        findex (findex < 1) = 1;                                           % Added for rebuttal


        % Smoothing of ridges
        beta1 = 1e-5; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
        iniIF = curvesmooth(Af(findex),beta1); % the initial guess for the IFs by ridge detection; the initial IFs should be smooth and thus we smooth the detected ridges

        % run the CT with g for TFC domain
        alpha1 = 50;
        tt = -8:1/Hz:8;
        h1 = exp(-pi*alpha1*tt.^2'); % window g_0
        eta = 1.5;
        [tfc1, tfrtic, tcrtic] = sqSTCTm(Sig, eta/length(Sig), 1, h1);

        % TFC ridges
        if1 = iniIF(1,:);
        if2 = iniIF(2,:);
        chirpn1 = gradient(if1)/(t(2)-t(1));
        chirpn2 = gradient(if2)/(t(2)-t(1));

        line1 = ones(length(Sig),Num); line2 = line1;                     
        for i = 1:length(t)
            [~,bb] = min(abs(tfrtic-if1(i)/Hz));
            line1(i,1) = bb;
            [~,bb1] = min(abs(tcrtic-chirpn1(i)/Hz^2));
            line1(i,2) = bb1-1;
            [~,bb2] = min(abs(tfrtic-if2(i)/Hz));
            line2(i,1) = bb2;
            [~,bb3] = min(abs(tcrtic-chirpn2(i)/Hz^2));
            line2(i,2) = bb3-1;
        end
        
        line1 (line1 < 1) = 1;                                                 % Added for rebuttal
        line2 (line2 < 1) = 1;                                                 % Added for rebuttal

        rrcon = zeros(2,length(Sig));
        scale = alpha1;
        for char = 1:length(Sig)
            tmp = zeros(2,2);
            aa=scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn1(char));
            tmp(1,1) = 1/sqrt(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn1(char)))*exp(-pi*(tfrtic(line1(char,1))*Hz-if1(char))^2/(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn1(char))));
            tmp(1,2) = 1/sqrt(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn2(char)))*exp(-pi*(tfrtic(line1(char,1))*Hz-if2(char))^2/(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn2(char))));
            tmp(2,1) = 1/sqrt(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn1(char)))*exp(-pi*(tfrtic(line2(char,1))*Hz-if1(char))^2/(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn1(char))));
            tmp(2,2) = 1/sqrt(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn2(char)))*exp(-pi*(tfrtic(line2(char,1))*Hz-if2(char))^2/(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn2(char))));

            xtmp = [tfc1(line1(char,2),line1(char,1),char); tfc1(line2(char,2),line2(char,1),char)]/Hz;
            rrcon(:,char) = tmp\xtmp;
        end

        % Rebuttal

        Sig1_s = Sig1(100: end-100); Sig2_s = Sig2(100: end-100);
        mode = rrcon(1,:)'; mode2 = rrcon(2,:)';
        mode_s = mode(100: end-100); mode2_s = mode2(100: end-100);


        e1 = norm(real(Sig1_s-mode_s))/norm(real(mode_s));
        e2 = norm(real(Sig2_s-mode2_s))/norm(real(mode2_s));

        e_k(k) = mean([e1 e2]);
    end
    e(j) = mean(e_k);
    j
end
figure;plot(real(Sig1_s));hold on;plot(real(mode_s),'r--')
figure;plot(real(Sig2_s));hold on;plot(real(mode2_s),'r--')
e=e';
% figure;plot(e)
toc

clear;clc;close all

Hz = 400; tend = 2;
t = 1/Hz : 1/Hz : tend; t = t'; Sf = 100;% time = time(1:400);
% 
amp1 = exp(0.1*t);  amp2 = exp(0.2*t);  amp3 = exp(0.3*t);
% 
Sig1 = amp1.*exp(2*1i*pi*(25*t.^2));                        IF1 = 50*t;
Sig2 = amp2.*exp(2*1i*pi*(25*t));                           IF2 = 25*ones(length(t),1);
Sig3 = amp3.*exp(2*1i*pi*(-22*t.^2 + 95*t + sin(5*t)));     IF3 = -44*t + 95 + 5*cos(5*t);

figure; set(gcf,'Position',[20 100 800 450]);
subplot(221);plot(t,amp1,'b','LineWidth',2);hold on;plot(t,amp2,'r','LineWidth',2);plot(t,amp3,'Color','[0, 0.5, 0]','LineWidth',2);
axis xy; xlim([0,t(end)]);ylim([0.9 1.9]);
legend('Mode 1','Mode 2','Mode 3','Location','northwest')
xlabel('Time (s)','FontSize',12,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',12,'FontName','Times New Roman');

j = 0; pin = 0; 
for SNR = -4:1:15
    j = j + 1;
    for k = 1 : 10
        pin = pin + 1;
        x = Sig1 + Sig2+ Sig3;    x = x';                                        % transposed
        seed = pin; % Use the loop variable as the seed
        rng(seed); % Set the seed for random number generation
        x = awgn(x,SNR);

        ng = 151;
        Nx = length(x);
        nw_2 = Nx;
        ntau = Nx;
        nc = 71;

        L = length(x);
        nw = round(nw_2/2);
        nc = nc+1-rem(nc,2);
        ng = ng+1-rem(ng,2);
        ng_h = (ng-1)/2;
        nc_2 = nc*2-1;
        cend = Hz*Hz/ng;

        % Atau, Aw and Ac are arguments of LCT time, frequency and chirp rate, respectively.
        Atau = (1/ntau:1/ntau:1)*L/Hz;
        Anw_2 = (0:(nw_2-1))*Hz/nw_2;
        Aw = Anw_2(1:nw);
        Ac = (-1:2/(nc-1):1)*cend;
        %
        % Calculation of signal LCT
        g = 0.54 - 0.46*cos(2.0*pi*(1:ng)'/(ng+1));  % HAMMING window
        g = g/sum(g);
        taon = round(Atau*Hz);
        x = hilbert(x);
        xext = [zeros(1,ng_h) x zeros(1,ng_h) ];     % extended signal
        xmat = xext(repmat((0:ng-1),[ntau,1,nc])+repmat(taon.',[1,ng,nc]));% matrix of truncated signal
        gmat = repmat(g.',[ntau,1,nc]);
        gt = (-ng_h:ng_h)/Hz;                        % time of window
        kn(1,:,:) = exp(-1j*pi*(gt).^2.'*Ac);        % kernel of LCT
        knmat = repmat(kn,[ntau,1,1]);
        xLCTa = xmat.*gmat.*knmat;
        xLCTb = zeros(ntau,nw_2,nc);
        xLCTb(:,1:ng,:) = xLCTa;
        xLCTc = circshift(xLCTb,-ng_h,2);
        xLCTc = fft(xLCTc,[],2);
        xLCT = xLCTc(:,1:nw,:);

        % Calculation of LCT Reference
        Ac_2 = (-2:2/(nc-1):2)*cend;

        RefLCTb = zeros(nw_2,nc_2);
        Refhmat = repmat(g,[1,nc_2]);
        Refknmat = exp(-1j*pi*(gt).^2.'*Ac_2);
        RefLCTa = Refhmat.*Refknmat;
        RefLCTb(1:ng,:) = RefLCTa;
        RefLCTc = circshift(RefLCTb,-ng_h,1);
        RefLCTc = fft(RefLCTc,[],1);
        LCTRef = fftshift(RefLCTc,1);

        fmidn = round((nw_2+1)/2)*ones(1,ntau); % frequency center index of LCT Reference
        cmidn = (nc_2+1)/2*ones(1,ntau);        % chirp rate center index of LCT Reference

        %iterative process
        TFR = zeros(ntau,nw);                       %% ---- for the other component, run this line again ---- %%%
        mf = 0 : nw-1;
        mc(1,1,:) = (0:nc-1)*nw_2;
        mfc = repmat(mf,ntau,1,nc)+repmat(mc,ntau,nw,1);
        % xenergy = sum(sum(xLCT(:,:,(nc+1)/2).*conj(xLCT(:,:,(nc+1)/2))));%signal energy
%         xenergy = sum(sum(xLCT(:,:,(nc+0)/1).*conj(xLCT(:,:,(nc+0)/1))));%signal energy

        % ratio = 1;
        % while ratio>threshold
        % Step 2
        absxLCT = abs(xLCT);
        [~,TFCMaxIndex] = max(absxLCT,[],[2,3],'linear');
        [~,maxfn,maxcn] = ind2sub([ntau,nw,nc],TFCMaxIndex);   % positions in matrix % CL (3,1);, (2,2);, ...
        TFMaxIndex = sub2ind([ntau,nw],1:ntau,maxfn.');
        mode = xLCT(TFCMaxIndex.');
        TFR(TFMaxIndex) = TFR(TFMaxIndex) + mode;

        % Step 3
        fdist = fmidn.'- maxfn + 1;                  % distance of frequency
        cdist = cmidn.'- maxcn + 1;                  % distance of chirp rate
        mtau = fdist + (cdist-1)*nw_2;
        m1b = repmat(mtau,1,nw,nc) + mfc;
        xLCTk = repmat(mode.',1,nw,nc).*LCTRef(m1b);  % LCT of the kth mode
        % Step 4
        xLCT = xLCT-xLCTk;% residual LCT
        % Step 5
        % renergy = sum(sum(xLCT(:,:,(nc+1)/2).*conj(xLCT(:,:,(nc+1)/2))));%residual energy
        % ratio = renergy/xenergy;
        % end
        % end

        TFR = zeros(ntau,nw);
        absxLCT = abs(xLCT);
        [~,TFCMaxIndex] = max(absxLCT,[],[2,3],'linear');
        [~,maxfn,maxcn] = ind2sub([ntau,nw,nc],TFCMaxIndex);   % positions in matrix % CL (3,1);, (2,2);, ...
        TFMaxIndex = sub2ind([ntau,nw],1:ntau,maxfn.');
        mode2 = xLCT(TFCMaxIndex.');
        TFR(TFMaxIndex) = TFR(TFMaxIndex) + mode2;
        % Step 3
        fdist = fmidn.'- maxfn + 1;                  % distance of frequency
        cdist = cmidn.'- maxcn + 1;                  % distance of chirp rate
        mtau = fdist + (cdist-1)*nw_2;
        m1b = repmat(mtau,1,nw,nc) + mfc;
        xLCTk = repmat(mode2.',1,nw,nc).*LCTRef(m1b);  % LCT of the kth mode
        % Step 4
        xLCT = xLCT-xLCTk;% residual LCT
        % Step 5
        % renergy = sum(sum(xLCT(:,:,(nc+1)/2).*conj(xLCT(:,:,(nc+1)/2))));%residual energy
        % ratio = renergy/xenergy;

        TFR = zeros(ntau,nw);
        absxLCT = abs(xLCT);
        [~,TFCMaxIndex] = max(absxLCT,[],[2,3],'linear');
        [~,maxfn,maxcn] = ind2sub([ntau,nw,nc],TFCMaxIndex);   % positions in matrix % CL (3,1);, (2,2);, ...
        TFMaxIndex = sub2ind([ntau,nw],1:ntau,maxfn.');
        mode3 = xLCT(TFCMaxIndex.');
        TFR(TFMaxIndex) = TFR(TFMaxIndex) + mode3;
        % Step 3
        fdist = fmidn.'- maxfn + 1;                  % distance of frequency
        cdist = cmidn.'- maxcn + 1;                  % distance of chirp rate
        mtau = fdist + (cdist-1)*nw_2;
        m1b = repmat(mtau,1,nw,nc) + mfc;
        xLCTk = repmat(mode3.',1,nw,nc).*LCTRef(m1b);  % LCT of the kth mode
        % Step 4
        xLCT = xLCT-xLCTk;% residual LCT
        %%%%%%%%%%   Rebuttal    %%%%%%%%%%
        % From VNCMD, eq. 40

        chops = 100;
        Sig1_s = Sig1(chops: end-chops); Sig2_s = Sig2(chops: end-chops); Sig3_s = Sig3(chops: end-chops);
        mode_s = mode(chops: end-chops); mode2_s = mode2(chops: end-chops); mode3_s = mode3(chops: end-chops);  


        e1 = norm(real(Sig1_s-mode_s'))/norm(real(mode_s'));
        e2 = norm(real(Sig2_s-mode2_s'))/norm(real(mode2_s'));
        e3 = norm(real(Sig3_s-mode3_s'))/norm(real(mode3_s'));

        e1_r = norm(real(Sig2_s-mode_s'))/norm(real(mode_s'));               % Because, may Sig1_s can be displayed as Sig2_s
        e2_r = norm(real(Sig3_s-mode2_s'))/norm(real(mode2_s'));
        e3_r = norm(real(Sig1_s-mode3_s'))/norm(real(mode3_s'));

        e1_r1 = norm(real(Sig3_s-mode_s'))/norm(real(mode_s'));               % Because, may Sig1_s can be displayed as Sig2_s
        e2_r1 = norm(real(Sig1_s-mode2_s'))/norm(real(mode2_s'));
        e3_r1 = norm(real(Sig2_s-mode3_s'))/norm(real(mode3_s'));
        
        e_M1 = [e1 e1_r e1_r1];
        e_M2 = [e2 e2_r e2_r1];
        e_M3 = [e3 e3_r e3_r1];

        e1_s = min(e_M1); e2_s = min(e_M2); e3_s = min(e_M3);

        e_k(k) = mean([e1_s e2_s e3_s]);
%         k
    end
    e(j) = mean(e_k);
    j
end
figure;
subplot(311);plot(real(Sig1_s));hold on;plot(real(mode_s),'r--')
subplot(312);plot(real(Sig2_s));hold on;plot(real(mode2_s),'r--')
subplot(313);plot(real(Sig3_s));hold on;plot(real(mode3_s),'r--')

% figure;plot(-25:1:15,e)
e = e';

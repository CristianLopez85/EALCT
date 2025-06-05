function [TFR,Atau,Af] = mALCT(x,fs,Nh,Nc,num)
%Adaptive linear chirplet transform (ALCT) will be soon published on IEEE transactions on industrial electronics
%Nx, Nf,Ntau, Nc, Hh are length of discrete signal, linear chirplet transform (LCT) frequency, LCT time, LCT chirp rate and window, respectively.
%Nfext and Ncext are about two times of NF and Nc.

Nx=length(x); Nfext=Nx; Ntau=Nx;

Nf=round(Nfext/2);
Nc=Nc+1-rem(Nc,2);
Nh=Nh+1-rem(Nh,2);
Nhhalfs=(Nh-1)/2;
Ncext=Nc*2-1;
cend=fs*fs/Nh;

%Atau, Af and Ac are arguments of LCT time, frequency and chirp rate, respectively.
Afext=(0:(Nfext-1))*fs/Nfext;
Af=Afext(1:Nf);
Atau=(1/Ntau:1/Ntau:1)*Nx/fs;
Ac=(-1:2/(Nc-1):1)*cend;
Acext=(-2:2/(Nc-1):2)*cend;
% Ac = linspace(-200,200,Nc);         %%%%%%% These 2 modified by CL  %%%%%%%
% Acext = linspace(-400,400,2*Nc-1);

% Calculation of signal LCT
TLCTb=zeros(Ntau,Nfext,Nc);
h=0.54 - 0.46*cos(2.0*pi*(1:Nh)'/(Nh+1));%HAMMING window
h=h/sum(h);
taon=round(Atau*fs);
x=hilbert(x);
xext=[zeros(1,Nhhalfs) x zeros(1,Nhhalfs) ];% extended signal
ht=(-Nhhalfs:Nhhalfs)/fs;% time of window
xmat=xext(repmat((0:Nh-1),[Ntau,1,Nc])+repmat(taon.',[1,Nh,Nc]));% matrix of truncated signal
hmat=repmat(h.',[Ntau,1,Nc]);
kn(1,:,:)=exp(-1j*pi*(ht).^2.'*Ac);%kernel of LCT
knmat=repmat(kn,[Ntau,1,1]);
TLCTa=xmat.*hmat.*knmat;
TLCTb(:,1:Nh,:)=TLCTa;
TLCTc=circshift(TLCTb,-Nhhalfs,2);
TLCTc=fft(TLCTc,[],2);
TLCT=TLCTc(:,1:Nf,:);

% Calculation of LCT Reference
RefLCTb=zeros(Nfext,Ncext);
Refhmat=repmat(h,[1,Ncext]);
Refknmat=exp(-1j*pi*(ht).^2.'*Acext);
RefLCTa=Refhmat.*Refknmat;
RefLCTb(1:Nh,:)=RefLCTa;
RefLCTc=circshift(RefLCTb,-Nhhalfs,1);
RefLCTc=fft(RefLCTc,[],1);
LCTRef=fftshift(RefLCTc,1);
fmidn=round((Nfext+1)/2)*ones(1,Ntau);% frequency center index of LCT Reference
cmidn=(Ncext+1)/2*ones(1,Ntau);% chirp rate center index of LCT Reference

%iterative process
TFR=zeros(Ntau,Nf);
mf=0:Nf-1;
mc(1,1,:)=(0:Nc-1)*Nfext;
mfc=repmat(mf,Ntau,1,Nc)+repmat(mc,Ntau,Nf,1);
% xenergy=sum(sum(xLCT(:,:,(Nc+1)/2).*conj(xLCT(:,:,(Nc+1)/2))));%signal energy
% ratio=1;
for i = 1 : num %while ratio>threshold
    absxLCT=abs(TLCT);
    [~,TFCMaxIndex]=max(absxLCT,[],[2,3],'linear');
    [~,maxfn,maxcn]=ind2sub([Ntau,Nf,Nc],TFCMaxIndex);
    mode=TLCT(TFCMaxIndex.');
    TFMaxIndex=sub2ind([Ntau,Nf],1:Ntau,maxfn.');
    TFR(TFMaxIndex)=TFR(TFMaxIndex)+mode;
    fdist=fmidn.'-maxfn+1;%distance of frequency
    cdist=cmidn.'-maxcn+1;%distance of chirp rate
    mtau=fdist+(cdist-1)*Nfext;
    m1b=repmat(mtau,1,Nf,Nc)+mfc;
    TLCTk=repmat(mode.',1,Nf,Nc).*LCTRef(m1b);% LCT of the kth mode
    TLCT=TLCT-TLCTk;% residual LCT

%     renergy=sum(sum(xLCT(:,:,(Nc+1)/2).*conj(xLCT(:,:,(Nc+1)/2))));%residual energy
    %     ratio=renergy/xenergy;
end
end

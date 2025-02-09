function [SNR]=calculate_SNR(AZT,DT,NA,Ntime,Stime)

band=0.2;
method=4;
j=1;i=1;
%############### SNR calculation ###############
% S波窗长
SWST=floor(Stime(i,1)/DT);
SWET=floor(Stime(i,2)/DT);
SDUL=SWET-SWST;
% Sdul=Stime(i,j+3)-Stime(i,j);
% 噪声窗长
NWST=floor(Ntime(i,1)/DT)+1;
NWET=floor(Ntime(i,2)/DT);
NDUL=NWET-NWST;
% Ndul=NWET*DT;
% 判断S 和 N 窗长是否需要分段计算SNR
if SDUL <= NDUL
    swave{1,j}(:,1)=AZT(SWST:SWET,j);
    noise{1,j}(:,1)=AZT(NWST:NWST+SDUL,j);
    swzct{1,j}(:,1)=ZC_TAPER(swave{1,j}(:,1),DT);
    nwzct{1,j}(:,1)=ZC_TAPER(noise{1,j}(:,1),DT);
    [~,~,freqs{1,j},SFASs{1,j}(:,1)]=FAS_Smooth(DT,swzct{1,j}(:,1),band,method);
    [~,~,freqs{1,j},NFASs{1,j}(:,1)]=FAS_Smooth(DT,nwzct{1,j}(:,1),band,method);
%     SNRraw{1,j}(:,1)=freq{1,j};
%     SNRraw{1,j}(:,2)=SFASraw{1,j}(:,1)./NFASraw{1,j}(:,1);
    SNR(:,1)=freqs{1,j};
    SNR(:,2)=SFASs{1,j}(:,1)./NFASs{1,j}(:,1);
else
    % S波分段窗
    kk1=fix(SDUL/NDUL);
    kk2=rem(SDUL,NDUL);
    % 整数部分以Noise长度为准
    noise{1,j}(:,1)=AZT(NWST:NWET,j);  % 整数部分
    nwzct{1,j}(:,1)=ZC_TAPER(noise{1,j}(:,1),DT);
    [freq{1,j},NFASraw{1,j}(:,1),freqs{1,j},NFASs{1,j}(:,1)]=FAS_Smooth(DT,nwzct{1,j}(:,1),band,method);
    for ss=1:kk1
        swaveseg{1,j}(:,ss)=AZT(SWST+(ss-1)*NDUL+1:SWST+ss*NDUL,j);
        swavesegzt{1,j}(:,ss)=ZC_TAPER(swaveseg{1,j}(:,ss),DT);
        [freq{1,j},swavesegFASraw{1,j}(:,ss),freqs{1,j},swavesegFAS{1,j}(:,ss)]=FAS_Smooth(DT,swavesegzt{1,j}(:,ss),band,method);
    end
    SsegFASraw{1,j}(:,1)=geomean(swavesegFASraw{1,j},2);
    SsegFAS{1,j}(:,1)=geomean(swavesegFAS{1,j},2);
    SNRraw1{1,j}(:,1)=SsegFASraw{1,j}(:,1)./NFASraw{1,j}(:,1); % 整数部分的SNR
    SNR1{1,j}(:,1)=SsegFAS{1,j}(:,1)./NFASs{1,j}(:,1);

    noiserem{1,j}(:,1)=AZT(NWST:NWST+kk2,j); % 余数部分
    nwremzct{1,j}(:,1)=ZC_TAPER(noiserem{1,j}(:,1),DT);
    nwremzcp{1,j}(:,1)=[nwremzct{1,j}(:,1);zeros(NDUL-kk2,1)];
    [freq{1,j},NRFASraw{1,j}(:,1),freqs{1,j},NRFASs{1,j}(:,1)]=FAS_Smooth(DT,nwremzcp{1,j}(:,1),band,method);
    swrem{1,j}(:,1)=AZT(SWET-kk2:SWET,j);
    swremzct{1,j}(:,1)=ZC_TAPER(swrem{1,j}(:,1),DT);
    swremzcp{1,j}(:,1)=[swremzct{1,j}(:,1);zeros(NDUL-kk2,1)];
    [freq{1,j},SRFASraw{1,j}(:,1),freqs{1,j},SRFAS{1,j}(:,1)]=FAS_Smooth(DT,swremzcp{1,j}(:,1),band,method);
    SNRraw1{1,j}(:,2)=SRFASraw{1,j}(:,1)./NRFASraw{1,j}(:,1); % 余数部分的SNR
    SNR1{1,j}(:,2)=SRFAS{1,j}(:,1)./NRFASs{1,j}(:,1);

    SNRraw{1,j}(:,1)=freq{1,j};
    SNR(:,1)=freqs{1,j};
%     SNRraw{1,j}(:,2)=geomean(SNRraw1{1,j},2);
    SNR(:,2)=geomean(SNR1{1,j},2);
end
end
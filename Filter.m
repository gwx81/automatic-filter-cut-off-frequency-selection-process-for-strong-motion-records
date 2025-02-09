function [ACC_sp,VEL_sp,DIS_sp,ACC,VEL,DIS,tend]=Filter(dt,acc,hpfc,lpfc,type,baseTF)

% padding
accraw(:,1) = acc(:,1)./980;
npts = length(accraw(:,1));
totallength = 2^(ceil(log2(npts))+1);
tbegin=ceil((totallength-npts)/2);
tend=totallength-tbegin-npts;
accpad(:,1)=[zeros(tbegin,1);accraw(:,1);zeros(tend,1)];
% FAS求解
FASf = fft(accpad)*dt;
df = 1/(totallength*dt);
Nnyq = totallength/2 + 1;
freq = ((1:Nnyq) - 1) * df;
Npts=length(freq);
FASFF=FASf(1:Npts,1);
FAS00 = abs(FASf);
FAS0(1:Npts,1) =FAS00(1:Npts,1);
% 滤波
nPole_HP=5;
nPole_LP=4;

% Acausal and causal Butterworth filter
if isnan(hpfc)==0 && hpfc>0
    butterworth = fbutter(Nnyq,df,hpfc,nPole_HP,'h',type);
    for kk=1:length(FASFF)
        accfs(kk,1)=FASFF(kk,1)*butterworth(1,kk);
    end
else
    for kk=1:length(FASFF)
        accfs(kk,1)=FASFF(kk,1);
    end
end

FASFF=accfs;
if isnan(lpfc)==0 && lpfc>0
    butterworth = fbutter(Nnyq,df,lpfc,nPole_LP,'l',type);
    for kk=1:length(FASFF)
        accfs(kk,1)=FASFF(kk,1)*butterworth(1,kk);
    end
else
    for kk=1:length(FASFF)
        accfs(kk,1)=FASFF(kk,1);
    end
end

accfs(1,1)= complex(real(accfs(1,1)),0);
accfs(end,1)= complex(real(accfs(end,1)),0);
accfs(totallength+2-(2:totallength/2))=conj(accfs(2:totallength/2));
ACC = ifft(accfs/dt)*980;
ACC(1,1)=0;

%baseline correction
if baseTF=='T'
    [ACC_sp,VEL_sp,DIS_sp]=BslnAdj(ACC,dt,tbegin,tend);
else
    
    % 积分
    VEL=cumtrapz(ACC)*dt;
    DIS=cumtrapz(VEL)*dt;
    ACC_sp=ACC(tbegin:length(ACC)-tend-1,1);
    VEL_sp=VEL(tbegin:length(ACC)-tend-1,1);
    DIS_sp=DIS(tbegin:length(ACC)-tend-1,1);
    % GRAVITY=980.665;
    % VEL=VEL*GRAVITY;
    % DIS=DIS*GRAVITY;
end



end
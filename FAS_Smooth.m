function [freq,FAS0,freq1,FAS]=FAS_Smooth(dt,acc,band,method)

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
FAS00 = abs(FASf);
FAS0(1:Npts,1) =FAS00(1:Npts,1);
% FAS = abs(FASf);
% N=length(FAS);
if band ==0
    return
end

if method == 1
    %FAS 平滑 parzen window
    if band==0
        return
    end
    T=1/df;
    UDF=1.854305/band*df;
    Lmax=floor(2/UDF)+1;
    if UDF > 0.5
        disp('BANDWIDTH IS TOO NARROW');
        while UDF <= 0.5
            band=band+0.05;
            UDF=1.854305/band*df;
            Lmax=floor(2/UDF)+1;
        end
        
    end
    if Lmax>101
        disp('BANDWIDTH IS TOO WIDE');
        while Lmax<=101
            band=band-0.05;
            UDF=1.854305/band*df;
            Lmax=floor(2/UDF)+1;
        end
    end
    
    % SPECTRAL WINDOW
    W(1,1)=0.75*UDF;
    for L=2:Lmax
        DIF=1.570796*single(L-1)*UDF;
        W(L,1)=W(1,1)*(sin(DIF)/DIF)^4;
    end
    % CONVERSION FROM FOURIER TO POWER SPECTRUM
    G(1,1)=FAS(1,1)^2/T;G(N,1)=FAS(N,1)^2/T;
    for K=2:N-1
        G(K,1)=2*FAS(K,1)^2/T;
    end
    % SMOOTHING OF POWER SPECTRUM
    if band~=0
        LL=Lmax*2-1;
        LN=LL-1+N;
        LT=(LL-1)*2+N;
        LE=LT-Lmax+1;
        G1=zeros(LT,1);
        G2=zeros(LE,1);
        for K=1:N
            G1(LL-1+K,1)=G(K,1);
        end
        for K=Lmax:LE
            S=W(1,1)*G1(K,1);
            for L=2:Lmax
                S=S+W(L,1)*(G1(K-L+1,1)+G1(K+L-1,1));
            end
            G2(K,1)=S;
        end
        for L=2:Lmax
            G2(LL+L-1,1)=G2(LL+L-1,1)+G2(LL-L+1,1);
            G2(LN-L+1,1)=G2(LN-L+1,1)+G2(LN+L-1,1);
        end
        for K=1:N
            G(K,1)=G2(LL-1+K,1);
        end
        %  SMOOTHED FOURIER SPECTRUM
        FAS(1,1)=sqrt(G(1,1)*T);
        FAS(N,1)=sqrt(G(N,1)*T);
        for K=2:N-1
            FAS(K,1)=sqrt(G(K,1)*T/2);
        end
        
    else
        return
    end
end

if method == 2
    %FAS 平滑 Konno window
    if band==0
        return
    end
    
    pi= 4.0*atan(1.0);
    work=zeros(Npts,1);
    for jj=1:Npts
        work(jj,1)=FAS(jj,1);
    end
    
    b=2.0*pi/band;
    %     b=4.0/band;
    FAS=zeros(Npts,1);
    f_fc=0.01^(1/b);
    for ic=1:Npts
        iz1 = round(ic*f_fc);
        iz2 = round(ic/f_fc);
        %         iz1 = round(ic*10.0^(-1.5*band));
        %         iz2 = round(ic*10.0^(1.5*band));
        %         % ractangle window
        %         iz1 = ic;
        %         iz2 = ic+10;
        if iz1 < 1
            iz1 = 1;
        end
        if iz2 > Npts
            iz2 = Npts;
        end
        sum = 0;
        wavg =0;
        for jj=iz1:iz2
            weight = w_konno_ohmachi(jj, ic, b);
            wavg = wavg + weight;
            sum = sum + weight*work(jj,1);
        end
        
        FAS(ic,1)=sum/wavg;
    end
    
end

if method == 3
    %FAS 平滑 五点三次平滑法
    mk = 15;
    for k=1:mk
        b(1,1)=(69*FAS(1,1)+4*(FAS(2,1)+FAS(4,1))-6*FAS(3,1)-FAS(5,1))/70;
        b(2,1)=(2*(FAS(1,1)+FAS(5,1))+27*FAS(2,1)+12*FAS(3,1)-8*FAS(4,1))/35;
        for jm=3:N-2
            b(jm,1)=(-3*(FAS(jm-2,1)+FAS(jm+2,1))+12*(FAS(jm-1,1)+FAS(jm+1,1))+17*FAS(jm,1))/35;
        end
        b(N-1,1)=(2*(FAS(N,1)+FAS(N-4,1))+27*FAS(N-1,1)+12*FAS(N-2,1)-8*FAS(N-3,1))/35;
        b(N,1)=(69*FAS(N,1)+4*(FAS(N-1,1)+FAS(N-3,1))-6*FAS(N-2,1)-FAS(N-4,1))/70;
        FAS(:,1)=b(:,1);
    end
    
end

if method == 4
    %FAS 平滑 矩形平滑法
    Y=FAS0(:,1);
    NN=Npts;
    NS=51;
    FS=freq';
    FMIN=min(FS);
    FMAX=max(FS);
    key1=-0.05;
    if key1>0
        W = zeros(abs(NS),1);
        S = zeros(NN,1);
        F = zeros(NN,1);
    end
    MAXNS=100;
    %     RAWVAR=0;
    %     RNDF=0;
    if key1<=0
        DF=FS(2,1)-FS(1,1);
        FINC=abs(key1);
        FSTART = log10(max(FMIN,0.001));
        [AVG,place]=SMUTH2(Y,NN,FMAX,FSTART,FINC,DF);
        
        [~,B]=findpeaks(place);
        if isempty(B)
            if place(end)>place(1)
                B=length(place);
            else
                B=1;
            end
        end
        if place(B,1)>length(freq)
            for jc=1:B-1
                freq1(jc,1)=FS(ceil(place(jc,1)),1);
                FAS(jc,1)=AVG(jc,1);
            end
            freq1(B,1)=FS(end,1);
            FAS(B,1)=FAS0(end,1);
        else
            for jc=1:B
                freq1(jc,1)=FS(ceil(place(jc,1)),1);
                FAS(jc,1)=AVG(jc,1);
            end
%             freq1(B+1,1)=FS(end,1);
%             FAS(B+1,1)=FAS0(end,1);
        end
        
        %         NFS=length(FAS);
        %         for ic=1:NFS
        %             F(ic,1)=place(ic);
        %         end
    else
        %基本上没用
        if NS>MAXNS
            print ('NS IS TOO LARGE, AUTO RESET TO ');
            print (MAXNS);
            NS=MAXNS;
        end
        if rem(NS,2)==0
            print ('EXITING IN SMOOTH: NS MUST BE ODD!!');
            return
        end
        NNS = (abs(NS)+1)/2;
        NP = NN+abs(NS)/2;
        NSAB = abs(NS);
        if NS>0
            % 矩形平滑
            for ic=1:NSAB
                W(ic,1)=1./NSAB;
            end
        else
            % 三角平滑
            ISUM =0;
            for ic= 1:NNS-1
                ISUM=ISUM+ic;
            end
            MS= NNS+2*ISUM;
            for ic=1:NNS
                J=NSAB-ic+1;
                W(ic,1)=ic/MS;
                W(J,1)=W(ic,1);
            end
        end
        % 平滑谱
        for ic=1:NN
            S(ic,1)=Y(ic,1)*W(NNS,1);
            for jc= 1:NNS-1
                if ic-jc<1
                    A=0;
                else
                    A=Y(ic-jc,1);
                end
                if ic-jc>NP
                    B=0;
                else
                    B=Y(ic+jc,1);
                end
                S(ic,1)=S(ic,1)+W(NNS-jc)*A+W(NNS+jc)*B;
            end
        end
        FAS=S;
    end
end

end
function [AVG,PLACE]=SMUTH2(FFT,NPTS,FMAX,FL1,DFL,DF)
% AVG=zeros(NPTS,1);
% PLACE=zeros(NPTS,1);
IAV=0;
FL=FL1-DFL;
FBU=10.^FL1;

while FBU<FMAX
    IAV=IAV+1;
    FL =FL+DFL;
    FBL=10.^(FL-0.5*DFL);
    FBU=10.^(FL+0.5*DFL);
    INL = floor(FBL/DF)+2;
    if (INL-2)*DF==FBL
        INL=INL-1;
    end
    if INL<2
        INL=2;
    end
    INU= floor(FBU/DF)+1;
    if INU>NPTS
        INU = NPTS;
    end
    if INU<INL
        AVG(IAV,1)=9.99;
    end
    AVG(IAV,1)=0;
    PLACE(IAV,1)=(INU-INL)/2+INL;
    for jc=INL:INU
        AVG(IAV,1)=log10(FFT(jc,1))+AVG(IAV,1);
    end
    
    AVG(IAV,1)=AVG(IAV,1)/(INU-INL+1);
end
k=0;
for ic=1:IAV
    if isfinite(AVG(ic,1))==1
        k=k+1;
        AVG(k,1)=10^AVG(ic,1);
        PLACE(k,1)=PLACE(ic,1);
    end
end
end
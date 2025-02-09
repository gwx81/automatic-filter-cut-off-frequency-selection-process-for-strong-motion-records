function fLP=autochoice_fLP(acc,DT,NA,Pt,St,Se)
%Select fLP using SNR and the default value corresponding to sampling frequency
%   acc: acceleration time series
%   DT: the sampling interval equal to the reciprocal of the sampling frequency(fc)
%   NA: the number of points
%   Pt: the arrival time of P wave
%   St: the arrival time of S wave
%   Se: the end time of S wave

% DC-MEAN
ADC=DC_mean(acc);
% Zerocrossing_Taper
AZT=ZC_TAPER(ADC,DT);
%SNR calculation
SNR=calculate_SNR(AZT,DT,NA,[0,0.95*Pt],[St,Se]);
if DT<=0.005
    lp_max=70;
elseif DT==0.01
    lp_max=40;
elseif DT==0.05
    lp_max=20;
else
    lp_max=45;
end
fLP=lp_max;
for i = 2:length(SNR)
    if SNR(i,1)>10 && SNR(i,1)<lp_max
        if SNR(i-1,2)>= 3&& SNR(i,2)<=3
            x=[SNR(i-1,1),SNR(i,1)];
            y=[SNR(i-1,2)-3,SNR(i,2)-3];
            kk=polyfit(x,y,1);
            fLP=roundn((-kk(2)/kk(1)),0);
            break
        end
    end
end
end


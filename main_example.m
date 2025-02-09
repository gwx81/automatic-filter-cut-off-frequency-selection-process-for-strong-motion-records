%%% main_example.m is an example where the autoselection functions applied are autochoice_fLP.m and autochoice_fHP.m
clear
close all
clc
data=load("HC.DEYA..HNE.D.INT-20230517_0000187.ACC.CV.dat");%record
DT=0.01;%DT is the sampling interval equal to the reciprocal of the sampling frequency(fc)
NA=length(data);%NA is the number of points
mag=4.3;%mag is the magnitude of the earthquake(Ms)
direction=1;%direction is the direction of the record, EW:1, NS:2, UD:3

Pt=117.6;%Pt is the arrival time of P wave
St=144.8;%St is the arrival time of S wave
Se=250;%Se is the end time of S wave
acc=data-mean(data(1:Pt/DT));%acc is a base-corrected acceleration time series
%select fHP and fLP
fLP=autochoice_fLP(acc,DT,NA,Pt,St,Se)%Select the low-pass filter cutoff frequency
fHP=autochoice_fHP(acc,DT,NA,mag,Pt,direction,fLP)%Select the high-pass filter cutoff frequency
%Time series of filtering results(Acausal filtering)
    %ACC, VEL, DIS: Time series with zero-padding
    %ACC_sp, VEL_sp, DIS_sp: Time series without zero-padding
    %tend: The number of points in zero-padding
[ACC_sp,VEL_sp,DIS_sp,ACC,VEL,DIS,tend]=Filter(DT,acc,fHP,fLP,'A','F');
%%%Note that the filtering results at this time cannot be used for final
%%%engineering applications and should be processed in post-processing
%%%steps, such as Taper and 6th order polynomial removal trends

t=DT*(1:NA);
figure;
subplot(3,1,1)
plot(t,ACC_sp)
subplot(3,1,2)
plot(t,VEL_sp)
subplot(3,1,3)
plot(t,DIS_sp)
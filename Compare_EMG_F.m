%Compares the Force data to the EMG signal graphically
clearvars
close all
clc

%loading required packages
pkg load io
pkg load control
pkg load signal

Enames={['iBiodexEMG_P001_T001_NAIM0.csv'] ['iBiodexEMG_P001_T001_NAIM30.csv'] ['iBiodexEMG_P001_T001_NAIM60.csv'] ['iBiodexEMG_P001_T001_NAIM90.csv']};
Fnames={['P001_T000_I0_StaticOptimization_force.sto'] ['P001_T000_I30_StaticOptimization_force.sto'] ['P001_T000_I60_StaticOptimization_force.sto'] ['P001_T000_I90_StaticOptimization_force.sto']};
titles={['IM0'] ['IM30'] ['IM60'] ['IM90']};
figure
hold on
for i=1:4
  E=importdata(char(Enames(i)),',');
  F=importdata(char(Fnames(i)),'\t');
  Emax=max(E(:,7));
  Fmax=max(F.data(:,5));
  Enorm=E(:,7)./Emax;
  Fnorm=F.data(:,5)./Fmax;
  Fsize=length(F.data);
  subplot(2,2,i)
  plot(E(1:Fsize,1),Enorm(1:Fsize),E(1:Fsize,1),Fnorm)
  title(char(titles(i)))
  legend('EMG','SO Force')
endfor
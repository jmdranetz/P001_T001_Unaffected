clearvars
close all
clc

%loading required packages
pkg load io
pkg load control
pkg load signal


%imports goniometer data
data=importdata('C:\Users\jo801202\Desktop\OctaveWorkspace\Pilot Data\P001_T001\unaffected\GonioTot.csv',',');
srg=data.data(3,1)-data.data(2,1); %finds sampling rates 
dsr=.01;

%Interpolates and crops goniometer data to find collection interval
starts=[0 0 0 0 0 0 0 0];
ends=[0 0 0 0 0 0 0 0];
for i=1:length(data.colheaders)/2;
  k=1;
  for j=1:length(data.data)%# of data points in a given trial
    if ~isnan(data.data(j,i*2))&&data.data(j,i*2)~=0
      k=j;
    endif
  endfor 
  t=srg*(k-1); %time duration of trial based on # of data points
  y=interp1(0:srg:t,data.data(1:k,i*2),linspace(0,t,t/(dsr))); %interpolates data
  b=y';
  for j=201:length(b);
    if (abs(b(j)-mean(b(50:200))) > 2) && (starts(i) == 0)
      starts(i)=j-50;
    endif
  endfor
  b=b(starts(i):(end-50));
  ends(i)=length(b);
endfor
subn=inputdlg('Enter Subject Number','Input Subject # (3 digits)');
fid=fopen(['Rajagopal_Biodex_interval_P001_T',char(subn),'.csv'], "w");
dlmwrite(fid,[starts;ends],',');
fclose(fid);
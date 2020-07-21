clearvars
close all
clc

%loading required packages
pkg load io
pkg load control
pkg load signal

%Select filter variables
defaultans = {'20', '5', '2','2000'};
x = inputdlg({'Enter the high pass filter cutoff value','Enter the low pass filter cutoff value','Enter the order of Butterworth Filter','Enter EMG system sampling frequency'}, 'Input', [1 30; 1 30; 1 30;1 30], defaultans); 
hpvar = str2double(x{1});
lpvar = str2double(x{2});
bworder = str2double(x{3});
samfreq = str2double(x{4});

%imports each files intervals
ends=importdata('Rajagopal_Biodex_interval_P001_T001.csv',',');
dsr=.01;

%files to input
filenames={['BiodexEMG_P001_T001_NAIM0.csv'] ['BiodexEMG_P001_T001_NAIM30.csv'] ['BiodexEMG_P001_T001_NAIM60.csv'] ['BiodexEMG_P001_T001_NAIM90.csv'] ['BiodexEMG_P001_T001_NAIMOpt.csv'] ['BiodexEMG_P001_T001_NAIK120.csv'] ['BiodexEMG_P001_T001_NAIK240.csv'] ['BiodexEMG_P001_T001_NAIK360.csv']};
fnames={['iBiodexEMG_P001_T001_NAIM0.csv'] ['iBiodexEMG_P001_T001_NAIM30.csv'] ['iBiodexEMG_P001_T001_NAIM60.csv'] ['iBiodexEMG_P001_T001_NAIM90.csv'] ['iBiodexEMG_P001_T001_NAIMOpt.csv'] ['iBiodexEMG_P001_T001_NAIK120.csv'] ['iBiodexEMG_P001_T001_NAIK240.csv'] ['iBiodexEMG_P001_T001_NAIK360.csv']};
for f=1:length(filenames)
  data=importdata(char(filenames(f)),',');
  sre=data.data(3,1)-data.data(2,1);
  
  %High pass filter
  [bb,aa] = butter(bworder, hpvar/(samfreq/2),'high');
  datafilt=filtfilt(bb,aa,data.data(:,2:7));
  %Rectifying EMG data
  datafilt=abs(datafilt);
  %Low pass filter
  [bb,aa] = butter(bworder, lpvar/(samfreq/2),'low');
  datafilt=filtfilt(bb,aa,datafilt);
  %interpolation
  idata=[];
  %figure
  %hold on
  for i=1:6
    idata(:,i)=interp1(data.data(:,1),datafilt(:,i),linspace(ends(1,f)./100,(ends(1,f)+ends(2,f)-1)./100,ends(2,f)));
  %  plot(linspace(ends(1,f)./100,(ends(1,f)+ends(2,f))./100,ends(2,f)),idata(:,1)) %plots it
  endfor
  %legend('LH','RF','VL','VM','MG','MH')
  %saving file
  fid=fopen(char(fnames(f)),"w");
  dlmwrite(fid,[(linspace(0,ends(2,f)./100-.01,ends(2,f)))',idata]);
  fclose(fid);
endfor
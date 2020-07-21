clearvars
close all
clc

%loading required packages
pkg load io
pkg load control
pkg load signal

%imports biodex data
fdata=importdata('C:\Users\jo801202\Desktop\OctaveWorkspace\Pilot Data\P001_T001\unaffected\BiodexTot.csv',',');
%imports goniometer data
data=importdata('C:\Users\jo801202\Desktop\OctaveWorkspace\Pilot Data\P001_T001\unaffected\GonioTot.csv',',');
sr=fdata.data(3,1)-fdata.data(2,1); %finds sampling rates 
srg=data.data(3,1)-data.data(2,1);
dsr=.01;

%Interpolates and crops goniometer data
B=[0];
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
%  [bb,aa] = butter(2,5/(100/2),'low'); %Low pass filter
%  b=filtfilt(bb,aa,b);
  if (size(B(:,1)) == 1) %recombines data into a matrix padded with NaN
    B = b;
  elseif (length(b) <  length(B(:,1)))
    B = [B [b; NaN(length(B(:,1))-length(b),1)]];
  elseif (length(b) == length(B(:,1)))
    B = [B b];
  elseif (length(b) >  length(B(:,1)))
    B = [[B; NaN(length(b)-length(B(:,1)),length(B(1,:)))] b];
  end
endfor

%calibrating goniometer angles
zer=mean(B(5:50,1));
B=B-zer;
deg=mean(B(600:650,4))./90;
B=B./(deg.*1.6115);
%convert to radians
B=-B.*pi()./180;
  
%plotting results
figure
hold on
for i=1:8
  plot([dsr:dsr:(dsr.*ends(i))]',B(1:ends(i),i))
endfor
legend('I0','I30','I60','I90','IOpt','K120','K240','K360')

%asking user for subject # are creates corresponding .mot file names
subn=inputdlg('Enter Subject Number','Input Subject # (3 digits)');
type={'I0','I30','I60','I90','IOpt','K120','K240','K360'};
for i=1:length(ends)
  fnames{i}=['Biodex_kinematics_P001_T',char(subn),'_',char(type{i}),'.mot'];
endfor

%saving kinematic data
for i=1:length(ends)
  fid=fopen(char(fnames{i}), "w");
  title=[char(fnames{i});'version=1';['nRows=',num2str(ends(i))];'nColumns=2';'inDegrees=no';'endheader'];
  head=["time\tknee_angle_l"];
  dlmwrite(fid,title,'\0');
  dlmwrite(fid,head,'\0');
  dlmwrite(fid,[[0:dsr:dsr.*(ends(i)-1)]',-B(1:ends(i),i)],'\t');
  fclose(fid);
endfor


%calculating reaction force
A=[0];
l=.406;
mg=22.24;
fixA=[0 pi/6 pi/3 pi/2 pi*49/180];
for i=1:8
  x=fdata.data(:,i*2); %data of a trial
  k=1;
  for j=1:length(x)%# of data points in a given trial
    if ~isnan(x(j))&&x(j)~=0
      k=j;
    endif
  endfor 
  t=sr*(k-1); %time duration of trial based on # of data points
  y=interp1(0:sr:t,x(1:k),linspace(0,t,t/(dsr))); %interpolates data
  b=y';
  b=b(starts(i):(end-50));
  if length(b) > length(B)
    b=b(1:length(B));
  endif
  Mb=6.0046.*b.^4+40.776.*b.^3+91.388.*b.^2-73.177.*b+9.5647; %finds moment about biodex arm
  if i <= 5
    Fr=Mb./l-mg.*sin((pi/2)-fixA(i))./2; %solves for reaction force isometric
  else
    Fr=Mb./l-mg.*sin((pi/2)-B(1:length(Mb),i))./2; %solves for reaction force
  endif
  if (size(A(:,1)) == 1) %recombines data into a matrix padded with NaN
    A = Fr;
  elseif (length(Fr) <  length(A(:,1)))
    A = [A [Fr; NaN(length(A(:,1))-length(Fr),1)]];
  elseif (length(Fr) == length(A(:,1)))
    A = [A Fr];
  elseif (length(Fr) >  length(A(:,1)))
    A = [[A; NaN(length(Fr)-length(A(:,1)),length(A(1,:)))] Fr];
  end
endfor

for i=1:8
  gnames{i}=['Biodex_reaction_force_P001_T',char(subn),'_',char(type{i}),'.mot'];
endfor

%saving reaction force data
for i=1:8
  fid=fopen(char(gnames{i}), "w");
  title=[char(gnames{i});'version=1';['nRows=',num2str(ends(i))];'nColumns=7';'inDegrees=no';'endheader'];
  head=["time\treaction_force_vx\treaction_force_vy\treaction_force_vz\treaction_force_px\treaction_force_py\treaction_force_pz"];
  dlmwrite(fid,title,'\0');
  dlmwrite(fid,head,'\0');
  dlmwrite(fid,[[0:dsr:dsr.*(ends(i)-1)]',A(1:ends(i),i),ones(ends(i),1)*0,ones(ends(i),1)*0,ones(ends(i),1)*0,-ones(ends(i),1)*0.15,ones(ends(i),1)*0],'\t');
  fclose(fid);
endfor

%writes a file storing the starting and end points
%fid=fopen(['Rajagopal_Biodex_interval_P001_T',char(subn),'.csv'], "w");
%dlmwrite(fid,[starts;ends],',');
%fclose(fid);
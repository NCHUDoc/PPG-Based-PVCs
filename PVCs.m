%%  20170529 - show all data by PlotATM
clc;
close all;
clear;
plotATM('212m');

%% 20170530 - show 8 seconds data
clc;
close all;
clear;
Name ='212m';
infoName = strcat(Name, '.info');
matName = strcat(Name, '.mat');

load(matName);
fid = fopen(infoName, 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[freqint] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = freqint(2);
fgetl(fid);
    for i = 1:size(val, 1)
      [row(i), signal(i), gain(i), base(i), units(i)]=strread(fgetl(fid),'%d%s%f%f%s','delimiter','\t');
    end
 
fclose(fid);
val(val==-32768) = NaN;


for i = 1:size(val, 1)
    val(i, :) = (val(i, :) - base(i)) / gain(i);
end

x = (1:size(val, 2)) * interval;

x=x';
val=val';
% samples = 100000;
% plot(x(1+samples:1000+samples), val(1+samples:1000+samples,1));
% hold on
% plot(x(1+samples:1000+samples), val(1+samples:1000+samples,6),'r');
plot(x(1:1000), val(1:1000,1));
hold on
plot(x(1:1000), val(1:1000,6),'r');
% plot(x, val(:,1));
% hold on
% plot(x, val(:,6),'r');
for i = 1:length(signal)
    labels{i} = strcat(signal{i}, ' (', units{i}, ')'); 
end
legend(labels{1,1},labels{1,6});

xlabel('Time (sec)');
grid on;
%% Interpolation to 500Hz
% https://www.google.com.tw/url?sa=t&rct=j&q=&esrc=s&source=web&cd=13&ved=0ahUKEwjnkLL_npfUAhVDoZQKHRsmAjcQFgiAATAM&url=https%3A%2F%2Fmirlab.org%2Fjang%2Fbooks%2FmatlabProgramming4guru%2Fslide%2F09-%25E5%2585%25A7%25E6%258F%2592%25E6%25B3%2595.ppt&usg=AFQjCNF1yMD4McCGM7iTxCQuqiNDQAZo5Q
% http://blog.sina.com.cn/s/blog_4c7482f101009vm2.html
% downsample/resample/decimate/interp/rat
clc;
close all;
clear;
Name ='212m';
infoName = strcat(Name, '.info');
matName = strcat(Name, '.mat');

% samples = 100000;
load(matName);
fid = fopen(infoName, 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[freqint] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = freqint(2);
fgetl(fid);
    for i = 1:size(val, 1)
      [row(i), signal(i), gain(i), base(i), units(i)]=strread(fgetl(fid),'%d%s%f%f%s','delimiter','\t');
    end
 
fclose(fid);
val(val==-32768) = NaN;


for i = 1:size(val, 1)
    val(i, :) = (val(i, :) - base(i)) / gain(i);
end

x = (1:size(val, 2)) * interval;
x=x';
val=val';
% for j = 1:size(x, 1)-2
%     val1(j) = ((val(j+2,2) -val(j+1,2)) / (x(j+2) - x(j+1)))* (x(j) - x(j+1)) + val(j+1,2);
% end
val1 = interp(val(:,1),4);
val6 = interp(val(:,6),4);
x1 = interp(x,4);
% https://cn.mathworks.com/help/signal/ref/interp.html
% t = 0:0.001:1;
% x = sin(2*pi*30*t) + sin(2*pi*60*t);
% y = interp(x,4);
% 
% subplot 211
% stem(0:30,x(1:31),'filled','markersize',3)
% grid on
% xlabel 'Sample number',ylabel Original
% subplot 212
% stem(0:120,y(1:121),'filled','markersize',3)
% grid on
% xlabel 'Sample number',ylabel Interpolated

% range = 750000;  % 750000 samples = 750000x0.008(1/125)seconds = 6000 seconds
range = 1000; 
plot(x(1:range), val(1:range,1));
hold on
plot(x(1:range), val(1:range,6),'r');
hold on
% 
% plot(x1(1+(4*samples):4*(range+samples)), val1(1+(4*samples):4*(range+samples)),'cx');
% plot(x1(1+(4*samples):4*(range+samples)), val6(1+(4*samples):4*(range+samples)),'y+');
plot(x1(1:4*range), val1(1:4*range),'cx');
plot(x1(1:4*range), val6(1:4*range),'y+');
hold on
for i = 1:length(signal)
    labels{i} = strcat(signal{i}, ' (', units{i}, ')'); 
end
labels{1,8} = 'Interpo ECG';
labels{1,9} = 'Interpo PPG';
legend(labels{1,1},labels{1,3},labels{1,8},labels{1,9} );

xlabel('Time (sec)');
grid on;

hold off
% figure(2)
% plot(x(1:1000), val1(1:1000),'c');
%% Preprocessing (filtering)
figure(1)
range = 6000*125*4; % interpolation (6000seconds(60min)*125(Hz)(samplerate)*4(after interpolation) =30000000)
% startpoint = 3600*125*4; % from 1hr start to 100 minutes(6000seconds) later = 3600seconds(60min)*125(Hz)(samplerate)*4(after interpolation) =  1800000
% val16000s = val1(startpoint:startpoint+range);
% val66000s = val6(startpoint:startpoint+range);

val16000s = val1(1:range);
val66000s = val6(1:range);
x16000s = x1(1:range);
plot(x16000s(1:1000),val16000s(1:1000));
hold on;
plot(x16000s(1:1000),val66000s(1:1000),'r');
grid on;
xlabel('Time (sec)');
% hold off;
% filter 


figure(2)

% Method 1:

% fs=500;		% Sampling rate
% filterOrder=2;		% Order of filter
% cutOffFreq1=0.001;	% Cutoff frequency
% cutOffFreq2=0.007;	% Cutoff frequency
% [b, a]=butter(filterOrder, [cutOffFreq1 cutOffFreq2], 'bandpass');
% % === Plot frequency response
% [h, w]=freqz(b, a);
% plot(w/pi*fs/2, abs(h), '.-'); title('Magnitude frequency response');
% grid on

% Method 2:
b = fir1(16,[0.00055 0.22223]);
freqz(b)

figure(3)
B= val16000s; 
val16000sf = filter(b,1,B);
Bp= val66000s; 
val66000sf = filter(b,1,Bp);
plot(x16000s(1:3000),val16000s(1:3000));
hold on;
plot(x16000s(1:3000),val66000s(1:3000),'r');
hold on;
plot(x16000s(1:3000),val16000sf(1:3000),'c');
hold on;
plot(x16000s(1:3000),val66000sf(1:3000),'y');
grid on;
xlabel('Time (sec)');
ylabel('Amp (mv)');

%% find beats

% Method 1:
Fs=500;
for i=1:30
    Mso_chan1(val16000sf(1+(i-1)*100000:i*100000),Fs) 
    ansx{i} = ans{:};
end

for i=1:30
beat(i) = length(ansx{i})
end
sum(beat)
% 
% so_chan(val16000sf(1000001:2000000),Fs) 
% 
% so_chan(val16000sf(2000001:3000000),Fs) 
% 
% Fs=500;
%  i=30
%     so_chan(val16000sf(1+(i-1)*100000:i*100000),Fs) 
%     ansy{i} = ans{:};
% Method 2:

Mso_chan2(val16000sf,Fs) 

 % Method 3: 
 
 [pks,locs]=findpeaks(val16000sf,'MinPeakHeight',0.3);
 
fprintf('Find Peaks = %d\n',length(pks));

 %%  20170603 - test show seperate data by PlotATM 
clc;
close all;
clear;
figure(1)
plotATM('212mp');

hold on;
plotATM('212me');
hold off
%%   
clc;
close all;
clear;
Name ='212me';
infoName = strcat(Name, '.info');
matName = strcat(Name, '.mat');

load(matName);
fid = fopen(infoName, 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[freqint] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = freqint(2);
fgetl(fid);
    for i = 1:size(val)
      [row(i), signal(i), gain(i), base(i), units(i)]=strread(fgetl(fid),'%d%s%f%f%s','delimiter','\t');
    end
 
fclose(fid);
val(val==-32768) = NaN;


% for i = 1:size(val, 1)
    vale(i, :) = (val(i, :) - base(i)) / gain(i);
% end


x = (1:length(vale)) * interval;

x=x';
vale=vale';
% plot(x(1:1000), val(1:1000,1));
% hold on
% plot(x(1:1000), val(1:1000,3),'r');
plot(x, vale);
hold on
    ecg = strcat(signal{i}, ' (', units{i}, ')'); 
% legend(labels{1});
labels{1} =ecg;
xlabel('Time (sec)');
grid on;



Name2 ='212mp';
infoName2 = strcat(Name2, '.info');
matName2 = strcat(Name2, '.mat');
load(matName2);
fid2 = fopen(infoName2, 'rt');
fgetl(fid2);
fgetl(fid2);
fgetl(fid2);
[freqint2] = sscanf(fgetl(fid2), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval2 = freqint2(2);
fgetl(fid2);
    for i = 1:size(val)
      [row(i), signal2(i), gain(i), base(i), units(i)]=strread(fgetl(fid),'%d%s%f%f%s','delimiter','\t');
    end
 
fclose(fid);
val(val==-32768) = NaN;

 valp(i, :) = (val(i, :) - base(i)) / gain(i);
% end


xp = (1:length(valp)) * interval;

xp=xp';
valp=valp';
plot(xp, valp,'r');
hold off;
for i = 1:length(signal2)
    ppg = strcat(signal2{i}, ' (', units{i}, ')'); 
end
labels{2} = ppg;
legend(labels{1} ,labels{2});

%%

clc;
close all;
clear;
Name ='212me1hr';
infoName = strcat(Name, '.info');
matName = strcat(Name, '.mat');

load(matName);
fid = fopen(infoName, 'rt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
[freqint] = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
interval = freqint(2);
fgetl(fid);
    for i = 1:size(val)
      [row(i), signal(i), gain(i), base(i), units(i)]=strread(fgetl(fid),'%d%s%f%f%s','delimiter','\t');
    end
 
fclose(fid);
val(val==-32768) = NaN;


% for i = 1:size(val, 1)
    vale(i, :) = (val(i, :) - base(i)) / gain(i);
% end


x = (1:length(vale)) * interval;

x=x';
vale=vale';
% plot(x(1:1000), val(1:1000,1));
% hold on
% plot(x(1:1000), val(1:1000,3),'r');
plot(x, vale);
hold on
    ecg = strcat(signal{i}, ' (', units{i}, ')'); 
% legend(labels{1});
labels{1} =ecg;
xlabel('Time (sec)');
grid on;
legend(labels{1});

%%

Fs=125;

    Mso_chan1(vale,Fs) 

%     Fs=250;
Mso_chan2(vale,Fs) 

% [pks,locs] = findpeaks(data)        % Find peaks and their indices
[pks,locs]=findpeaks(vale,'MinPeakHeight',0.3);
  plot(vale,'Color','blue'); hold on;
%   plot(locs,vale(locs),'k^','markerfacecolor',[1 0 0]);
% ======================================================================================
%% WFDB tools from MIT - 2017060601  
% https://physionet.org/physiotools/matlab/wfdb-app-matlab/
%
% rdsamp read signal files of WFDB records
% ======================================================================================

clc 
clear
close all;
[tm, signal]=rdsamp('D:\MIT-BIH\MIT-BIH(Arrhythmia Database)\101',[],1000);
plot(tm,signal(:,1));
hold on;
plot(tm,signal(:,2),'r');
grid on;
%% -- 20170614 PM 1:42 --  Read data from MIT-BIH test%%
clc;
clear all;
%------ SPECIFY DATA -------%%%%%%  SPECIFY DATA的所有變數可視情況作調整-----
PATH= 'D:\MIT-BIH\MIT-BIH(Arrhythmia Database)';   % path, where data are saved
HEADERFILE= '101.hea';      % header-file in text format
ATRFILE= '101.atr';         % attributes-file in binary format
DATAFILE= '101.dat';        % data-file
% SAMPLES2READ=650000;        % number of samples to be read
% SAMPLES2READ=324000;
SAMPLES2READ=108000;
channel=1;                  % LeadII一般在1，Record114在2，Record102＆104僅有V5＆V2無LeadII。
% in case of more than one signal:
% 2*SAMPLES2READ samples are read
sample_rate=360;
%------ LOAD HEADER DATA --------------------------------------------------
fprintf(1,'WORKING ON %s ...\n', HEADERFILE);
signalh= fullfile(PATH, HEADERFILE);
fid1=fopen(signalh,'r');
z= fgetl(fid1);
A= sscanf(z, '%*s %d %d %d',[1,3]);
nosig= A(1);                % number of signals
sfreq= A(2);                % sample rate of data
clear A;
for kk=1:nosig
z= fgetl(fid1);
A= sscanf(z, '%*s %d %d %d %d %d',[1,5]);
dformat(kk)= A(1);           % format; here only 212 is allowed
gain(kk)= A(2);              % number of integers per mV
bitres(kk)= A(3);            % bitresolution
zerovalue(kk)= A(4);         % integer value of ECG zero point
firstvalue(kk)= A(5);        % first integer value of signal (to test for errors)
end;
fclose(fid1);
clear A;
%------ LOAD BINARY DATA --------------------------------------------------
if dformat~= [212,212], error('this script does not apply binary formats different to 212.'); end;
signald= fullfile(PATH, DATAFILE);            % data in format 212
fid2=fopen(signald,'r');
A= fread(fid2, [3, SAMPLES2READ], 'uint8')';  % matrix with 3 rows, each 8 bits long, = 2*12bit
fclose(fid2);
M2H= bitshift(A(:,2), -4);            % 字元向右移四位，即取字元的高四位
M1H= bitand(A(:,2), 15);              % 取字元的低四位
PRL=bitshift(bitand(A(:,2),8),9);     % sign-bit  取出字元低四位中最高位，向左移九位
PRR=bitshift(bitand(A(:,2),128),5);   % sign-bit  取出字元高四位中最高位，向左移五位
M( : , 1)= bitshift(M1H,8)+ A(:,1)-PRL;
M( : , 2)= bitshift(M2H,8)+ A(:,3)-PRR;
if M(1,:) ~= firstvalue, error('inconsistency in the first bit values'); end;
clear A M1H M2H PRR PRL;
fprintf(1,'LOADING DATA FINISHED \n');

ecgdata=M(:,channel);
ecg_data=ecgdata;

%% 
cmap = {'b', 'g','c','y','g'};
for channel = 1: 2
%     plot(M(:,channel),cmap{channel}); hold on
    plot(M(1:1000,channel),cmap{channel}); hold on
end
%% Read atr
[tm, signal]=rdsamp('D:\MIT-BIH\MIT-BIH(Arrhythmia Database)\107',[],650000);
[ann]=rdann('D:\MIT-BIH\MIT-BIH(Arrhythmia Database)\107', 'atr', [],[],[],'V');
figure(1)
plot(tm,signal(:,1));hold on;grid on
plot(tm(ann(:,1)),signal(ann(:,1),1),'ro','MarkerSize',8);
title('MIT-BIH(Arrhythmia Database) 107.dat')
xlabel('Samples');
ylabel('Amplitute');
legend('ECG','PVCs');
%% Count peak with MIT-BIH Arrhythmia DB signal 100 (5000 samples)
clear
close all
clc
[tm, signal]=rdsamp('D:\MIT-BIH\MIT-BIH(Arrhythmia Database)\100',[],5000);
plot(tm,signal(:,1));grid on

for i=1:length(signal(:,1))
    signalecg=signal(:,1);
end

indx=find(signalecg>0);

diffindx = indx(2:end) - indx(1:end -1); 
indgap=find(diffindx>1);

indmax=[]; % the location of index with maximal value in each cycle
for k=1:length(indgap)+1
    if k==1
        period=indx(1:indgap(1));
    elseif k==length(indgap)+1
        period=indx(indgap(k-1)+1:end);
    else
        period=indx(indgap(k-1)+1:indgap(k));
    end
    [value,ind]=max(signalecg(period));
    indmax(k)=period(ind(1));
end
figure,
plot(tm,signalecg);hold on;grid on
plot(indmax/360,signalecg(indmax), 'ro') 
title('MIT-BIH(Arrhythmia Database) 100.dat')
xlabel('Time(Seconds)');
ylabel('Amplitute');
legend('ECG','Peaks');


%% Read data from mimicdb
clc 
clear
close all;
[tm2, signal2]=rdsamp('D:\MIT-BIH\mimicdb\484\48400001',[],75000); % mimicdb maximum database 600 seconds
% plot(tm2,signal2(:,1));
% hold on;
% plot(tm2,signal2(:,7),'r');
plot(tm2(1:1000),signal2(1:1000,1));
hold on;
plot(tm2(1:1000),signal2(1:1000,7),'r');
%% 
clc 
clear
close all;
[tm3,signal3]=rdsamp('D:\MIT-BIH\mimicdb\482\48200003',[3 7],75000);
% plot(tm2,signal2(:,1));
% hold on;
% plot(tm2,signal2(:,7),'r');
% plot(tm3(8001:9000),signal3(8001:9000,1));
% hold on;
% plot(tm3(8001:9000),signal3(8001:9000,2),'r');
plot(tm3(4001:5000),signal3(4001:5000,1));
hold on;
plot(tm3(4001:5000),signal3(4001:5000,2),'r');
%% Test Table from reference paper for 9 signals with 1 to 2hr40min ECG and PPG datas
plotATM('482m');
%%
% Take 1-2hr PPG of signal 484, but ignore the gain value
clc 
clear
close all;

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\484\48400007',[2 7],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\484\48400008',[2 7],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\484\48400009',[2 7],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\484\48400010',[2 7],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\484\48400011',[2 7],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\484\48400012',[2 7],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\484\48400013',[2 7],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\484\48400014',[2 7],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\484\48400015',[2 7],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\484\48400016',[2 7],75000);

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\039\03900007',[1 3],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\039\03900008',[1 3],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\039\03900009',[1 3],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\039\03900010',[1 3],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\039\03900011',[1 3],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\039\03900012',[1 3],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\039\03900013',[1 3],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\039\03900014',[1 3],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\039\03900015',[1 3],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\039\03900016',[1 3],75000);

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\221\22100007',[1 4],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\221\22100008',[1 4],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\221\22100009',[1 4],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\221\22100010',[1 4],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\221\22100011',[1 4],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\221\22100012',[1 4],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\221\22100013',[1 4],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\221\22100014',[1 4],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\221\22100015',[1 4],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\221\22100016',[1 4],75000);

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\230\23000007',[1 5],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\230\23000008',[1 5],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\230\23000009',[1 5],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\230\23000010',[1 5],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\230\23000011',[1 5],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\230\23000012',[1 5],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\230\23000013',[1 5],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\230\23000014',[1 5],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\230\23000015',[1 5],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\230\23000016',[1 5],75000);

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\253\25300007',[2 6],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\253\25300008',[2 6],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\253\25300009',[2 6],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\253\25300010',[2 6],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\253\25300011',[2 6],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\253\25300012',[2 6],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\253\25300013',[2 6],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\253\25300014',[2 6],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\253\25300015',[2 6],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\253\25300016',[2 6],75000);

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\439\43900007',[1 6],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\439\43900008',[1 6],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\439\43900009',[1 6],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\439\43900010',[1 6],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\439\43900011',[1 6],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\439\43900012',[1 6],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\439\43900013',[1 6],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\439\43900014',[1 6],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\439\43900015',[1 6],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\439\43900016',[1 6],75000);

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\444\44400007',[1 5],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\444\44400008',[1 5],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\444\44400009',[1 5],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\444\44400010',[1 5],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\444\44400011',[1 5],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\444\44400012',[1 5],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\444\44400013',[1 5],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\444\44400014',[1 5],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\444\44400015',[1 5],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\444\44400016',[1 5],75000);

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\449\44900007',[1 5],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\449\44900008',[1 5],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\449\44900009',[1 5],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\449\44900010',[1 5],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\449\44900011',[1 5],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\449\44900012',[1 5],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\449\44900013',[1 5],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\449\44900014',[1 5],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\449\44900015',[1 5],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\449\44900016',[1 5],75000);

[tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\482\48200007',[3 7],75000);
[tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\482\48200008',[3 7],75000);
[tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\482\48200009',[3 7],75000);
[tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\482\48200010',[3 7],75000);
[tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\482\48200011',[3 7],75000);
[tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\482\48200012',[3 7],75000);
[tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\482\48200013',[3 7],75000);
[tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\482\48200014',[3 7],75000);
[tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\482\48200015',[3 7],75000);
[tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\482\48200016',[3 7],75000);

signal7ecg = signal7(:,1);
signal8ecg = signal8(:,1);
signal9ecg = signal9(:,1);
signalaecg = signala(:,1);
signalbecg = signalb(:,1);
signalcecg = signalc(:,1);
signaldecg = signald(:,1);
signaleecg = signale(:,1);
signalfecg = signalf(:,1);
signal1ecg = signal1(:,1);

signalecg484_1to2hr = [ signal7ecg' signal8ecg' signal9ecg' signalaecg' signalbecg' signalcecg' signaldecg' signaleecg' signalfecg' signal1ecg'];
signalecg484_1to2hr = signalecg484_1to2hr';

signal7ppg = signal7(:,2);
signal8ppg = signal8(:,2);
signal9ppg = signal9(:,2);
signalappg = signala(:,2);
signalbppg = signalb(:,2);
signalcppg = signalc(:,2);
signaldppg = signald(:,2);
signaleppg = signale(:,2);
signalfppg = signalf(:,2);
signal1ppg = signal1(:,2);

signalppg484_1to2hr = [ signal7ppg' signal8ppg' signal9ppg' signalappg' signalbppg' signalcppg' signaldppg' signaleppg' signalfppg' signal1ppg'];
signalppg484_1to2hr = signalppg484_1to2hr';

plot(signalecg484_1to2hr(1:1000));
hold on
plot(signalppg484_1to2hr(1:1000),'r');
%% Interpolation

signalecg484_1to2hrI = interp(signalecg484_1to2hr,4);
signalppg484_1to2hrI = interp(signalppg484_1to2hr,4);

%% Preprocessing (filtering)
b = fir1(16,[0.00055 0.22223]);
B= signalecg484_1to2hrI; 
signalecg484_1to2hrIf = filter(b,1,B);
B= signalppg484_1to2hrI; 
signalppg484_1to2hrIf = filter(b,1,B);

%% find beats

% Method 1:
Fs=500;
for i=1:30
    Mso_chan1(signalecg484_1to2hrIf(1+(i-1)*100000:i*100000),Fs) 
    ansx{i} = ans{:};
end


for i=1:30
beat(i) = length(ansx{i})
end
sum(beat)
% 
% so_chan(val16000sf(1000001:2000000),Fs) 
% 
% so_chan(val16000sf(2000001:3000000),Fs) 
% 
% Fs=500;
%  i=30
%     so_chan(val16000sf(1+(i-1)*100000:i*100000),Fs) 
%     ansy{i} = ans{:};
% Method 2:

Mso_chan2(signalecg484_1to2hrIf,Fs) 

 % Method 3: 
 
 [pks,locs]=findpeaks(signalecg484_1to2hrIf,'MinPeakHeight',0.3);
 
fprintf('Find Peaks = %d\n',length(pks));
%% Tranfer data for SingGUI use 20160612

clc 
clear
close all;
plotATM('484m')
close all
[tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\484\48400007',[1 7],75000);
[tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\484\48400008',[1 7],75000);
[tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\484\48400009',[1 7],75000);
[tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\484\48400010',[1 7],75000);
[tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\484\48400011',[1 7],75000);
[tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\484\48400012',[1 7],75000);
[tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\484\48400013',[1 7],75000);
[tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\484\48400014',[1 7],75000);
[tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\484\48400015',[1 7],75000);
[tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\484\48400016',[1 7],75000);

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\039\03900007',[1 3],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\039\03900008',[1 3],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\039\03900009',[1 3],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\039\03900010',[1 3],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\039\03900011',[1 3],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\039\03900012',[1 3],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\039\03900013',[1 3],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\039\03900014',[1 3],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\039\03900015',[1 3],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\039\03900016',[1 3],75000);

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\221\22100007',[1 4],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\221\22100008',[1 4],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\221\22100009',[1 4],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\221\22100010',[1 4],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\221\22100011',[1 4],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\221\22100012',[1 4],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\221\22100013',[1 4],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\221\22100014',[1 4],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\221\22100015',[1 4],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\221\22100016',[1 4],75000);

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\230\23000007',[1 5],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\230\23000008',[1 5],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\230\23000009',[1 5],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\230\23000010',[1 5],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\230\23000011',[1 5],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\230\23000012',[1 5],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\230\23000013',[1 5],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\230\23000014',[1 5],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\230\23000015',[1 5],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\230\23000016',[1 5],75000);

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\253\25300007',[2 6],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\253\25300008',[2 6],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\253\25300009',[2 6],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\253\25300010',[2 6],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\253\25300011',[2 6],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\253\25300012',[2 6],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\253\25300013',[2 6],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\253\25300014',[2 6],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\253\25300015',[2 6],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\253\25300016',[2 6],75000);

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\439\43900007',[1 6],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\439\43900008',[1 6],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\439\43900009',[1 6],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\439\43900010',[1 6],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\439\43900011',[1 6],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\439\43900012',[1 6],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\439\43900013',[1 6],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\439\43900014',[1 6],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\439\43900015',[1 6],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\439\43900016',[1 6],75000);

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\444\44400007',[1 5],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\444\44400008',[1 5],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\444\44400009',[1 5],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\444\44400010',[1 5],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\444\44400011',[1 5],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\444\44400012',[1 5],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\444\44400013',[1 5],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\444\44400014',[1 5],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\444\44400015',[1 5],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\444\44400016',[1 5],75000);

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\449\44900007',[1 5],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\449\44900008',[1 5],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\449\44900009',[1 5],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\449\44900010',[1 5],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\449\44900011',[1 5],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\449\44900012',[1 5],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\449\44900013',[1 5],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\449\44900014',[1 5],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\449\44900015',[1 5],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\449\44900016',[1 5],75000);

% [tm7,signal7]=rdsamp('D:\MIT-BIH\mimicdb\482\48200007',[3 7],75000);
% [tm8,signal8]=rdsamp('D:\MIT-BIH\mimicdb\482\48200008',[3 7],75000);
% [tm9,signal9]=rdsamp('D:\MIT-BIH\mimicdb\482\48200009',[3 7],75000);
% [tma,signala]=rdsamp('D:\MIT-BIH\mimicdb\482\48200010',[3 7],75000);
% [tmb,signalb]=rdsamp('D:\MIT-BIH\mimicdb\482\48200011',[3 7],75000);
% [tmc,signalc]=rdsamp('D:\MIT-BIH\mimicdb\482\48200012',[3 7],75000);
% [tmd,signald]=rdsamp('D:\MIT-BIH\mimicdb\482\48200013',[3 7],75000);
% [tme,signale]=rdsamp('D:\MIT-BIH\mimicdb\482\48200014',[3 7],75000);
% [tmf,signalf]=rdsamp('D:\MIT-BIH\mimicdb\482\48200015',[3 7],75000);
% [tm1,signal1]=rdsamp('D:\MIT-BIH\mimicdb\482\48200016',[3 7],75000);

signal7ecg = signal7(:,1);
signal8ecg = signal8(:,1);
signal9ecg = signal9(:,1);
signalaecg = signala(:,1);
signalbecg = signalb(:,1);
signalcecg = signalc(:,1);
signaldecg = signald(:,1);
signaleecg = signale(:,1);
signalfecg = signalf(:,1);
signal1ecg = signal1(:,1);

signalecg484_1to2hr = [ signal7ecg' signal8ecg' signal9ecg' signalaecg' signalbecg' signalcecg' signaldecg' signaleecg' signalfecg' signal1ecg'];
signalecg484_1to2hr = signalecg484_1to2hr';

signal7ppg = signal7(:,2);
signal8ppg = signal8(:,2);
signal9ppg = signal9(:,2);
signalappg = signala(:,2);
signalbppg = signalb(:,2);
signalcppg = signalc(:,2);
signaldppg = signald(:,2);
signaleppg = signale(:,2);
signalfppg = signalf(:,2);
signal1ppg = signal1(:,2);

signalppg484_1to2hr = [ signal7ppg' signal8ppg' signal9ppg' signalappg' signalbppg' signalcppg' signaldppg' signaleppg' signalfppg' signal1ppg'];
signalppg484_1to2hr = signalppg484_1to2hr';

plot(signalecg484_1to2hr(1:1000));
hold on
plot(signalppg484_1to2hr(1:1000),'r');
%%
close all
plot(signalecg484_1to2hr); grid on ; hold on

signalppg484_1to2hr(1:length(signalppg484_1to2hr)) = signalppg484_1to2hr(1:length(signalppg484_1to2hr)) -0.8;

plot(signalppg484_1to2hr,'r');
%%

ecg_1to2hr = [ signal7ecg signal8ecg signal9ecg signalaecg signalbecg signalcecg signaldecg signaleecg signalfecg signal1ecg];
ppg_1to2hr = [ signal7ppg signal8ppg signal9ppg signalappg signalbppg signalcppg signaldppg signaleppg signalfppg signal1ppg];
clear signal*

dlmwrite('484ecg.dat', ecg_1to2hr);
dlmwrite('484ppg.dat', ppg_1to2hr);

save('484ecg.mat', 'ecg_1to2hr');
save('484ppg.mat', 'ppg_1to2hr');

%% Test PPG peaks - 2017061401

clc
clear 
close all
% load ('484ecg.mat')
load ('484ppg.mat')
for i = 1:10
%     ecg{i}(1,:)= ecg_1to2hr(:,i);
   ecg{i}(1,:)= ppg_1to2hr(:,i);
end

% ecgtotal = [ecg{1}' ecg{2}' ecg{3}' ecg{4}' ecg{5}' ecg{6}' ecg{7}' ecg{8}' ecg{9}' ecg{10}'];
ecgtotal = [ecg{1} ecg{2} ecg{3} ecg{4} ecg{5} ecg{6} ecg{7} ecg{8} ecg{9} ecg{10}];
testtotal = ecgtotal';
 plot(testtotal,'c');
plot(ecgtotal(1:3000));

for j=0:74
    for i=1:length(ecgtotal(1,(1+j*10000):(10000+j*10000)))
        signalecg=ecgtotal(1,(1+j*10000):(10000+j*10000));
    end
    signalecg = signalecg';
    
    indx=find(signalecg>0.2);
    
    diffindx = indx(2:end) - indx(1:end -1);
    indgap=find(diffindx>1);
    
    indmax=[]; % the location of index with maximal value in each cycle
    for k=1:length(indgap)+1
        if k==1
            period=indx(1:indgap(1));
        elseif k==length(indgap)+1
            period=indx(indgap(k-1)+1:end);
        else
            period=indx(indgap(k-1)+1:indgap(k));
        end
        [value,ind]=max(signalecg(period));
        indmax(k)=period(ind(1));
    end
    peaks(j+1)=length(indmax);
end
figure,
plot(signalecg);hold on;grid on
plot(indmax,signalecg(indmax), 'ro') 
% title('MIT-BIH(Arrhythmia Database) 100.dat')
title('MIMIC DB Record 484')
xlabel('Time(Seconds)');
ylabel('Amplitute');
% legend('ECG','Peaks');
legend('PPG','Peaks');
peakstotal=sum(peaks)

%% Read MIMIC database annotation - 2017061601
clc
clear
close all
channel =1;
PATH = 'D:\MIT-BIH\mimicdb\484\';
DATAFILE = '484.qrs';
%------ LOAD BINARY DATA --------------------------------------------------
signald= fullfile(PATH, DATAFILE);            % data in format 212
fid2=fopen(signald,'r');
A= fread(fid2, 'uint8')';  % matrix with 3 rows, each 8 bits long, = 2*12bit
fclose(fid2);

fprintf(1,'LOADING DATA FINISHED \n');

char(A(5:16))

char(A(21:32))

% for i = 1:length(A)-15
%     QRS{i}=char(A(i*16-11:i*16))
% end
j=1;
for i = 1:length(A)-12
    if strcmp('Q',char(A(i)))
        if strcmp('R',char(A(i+1)))
            QRS{j}=char(A(i:i+12));
             j=j+1;
        end
     end
end
%%
% Subject 484 (m, 60) 
% Clinical class: 
% Record 484 Duration 44.9 hours 
% Signals: II I ABP PAP CVP LAP PLETH CO2 
% Measurements: ABP AWRR C.O. CVP ETCO2 HR IMCO2 LAP PAP SpO2 TBLOOD 
% HR_bpm = 88.303637713437270 average
time_hour=44.9;
time_min=44.9*60;
time_sec=44.9*60*60;
peaks=length(QRS);
HR_sec=peaks/time_sec;
HR_bpm=HR_sec*60 
HR_bph=HR_bpm*60 



%%
clear Str
for i=1:length(QRS)
    Str(i,1:3)=QRS{i}(1,6:8);
end
for i=1:length(Str)
    QRSw(i)=str2num(Str(i,1:3));
end
%%

bar(QRSw)
grid on
title('MIMIC database 484.qrs')
xlabel('Number of QRS (peaks)')
ylabel('QRS width (ms)')
legend('Annatation QRSw(ms)')

%% Get QRSw annotation from 1hour to 2hour 40 mins and analysis
% HR_bph = 4997

QRSw_1hr2hr40=QRSw(4997:(4997*2));
QRSw_avg=mean(QRSw_1hr2hr40)
abnormal =0;
abnormal2 =0;
j=1;k=1;
for i=1: length(QRSw_1hr2hr40)
    if (QRSw_1hr2hr40(i) - QRSw_avg)>5
        abnormal = abnormal +1;
        indx(j) = i;
        j=j+1;
    elseif(QRSw_1hr2hr40(i) - QRSw_avg)<=-6
        abnormal2 = abnormal2 +1;
        indx2(k) = i;
        k=k+1;
    end

end
%%



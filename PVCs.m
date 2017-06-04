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
plot(x(1:1000), val(1:1000,1));
hold on
plot(x(1:1000), val(1:1000,6),'r');

for i = 1:length(signal)
    labels{i} = strcat(signal{i}, ' (', units{i}, ')'); 
end
legend(labels{1,1},labels{1,6});

xlabel('Time (sec)');
grid on;
%% Interpolation to 500Hz
% https://www.google.com.tw/url?sa=t&rct=j&q=&esrc=s&source=web&cd=13&ved=0ahUKEwjnkLL_npfUAhVDoZQKHRsmAjcQFgiAATAM&url=https%3A%2F%2Fmirlab.org%2Fjang%2Fbooks%2FmatlabProgramming4guru%2Fslide%2F09-%25E5%2585%25A7%25E6%258F%2592%25E6%25B3%2595.ppt&usg=AFQjCNF1yMD4McCGM7iTxCQuqiNDQAZo5Q
% http://blog.sina.com.cn/s/blog_4c7482f101009vm2.html

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
startpoint = 3600*125*4; % from 1hr start to 100 minutes(6000seconds) later = 3600seconds(60min)*125(Hz)(samplerate)*4(after interpolation) =  1800000
val16000s = val1(startpoint:startpoint+range);
val66000s = val6(startpoint:startpoint+range);
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
 %%  20170603 - show seperate data by PlotATM
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

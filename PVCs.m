%%  20170529 - show all data by PlotATM
clc;
close all;
clear;
plotATM('039m');

%% 20170530 - show 8 seconds data
clc;
close all;
clear;
Name ='039m';
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
plot(x(1:1000), val(1:1000,3),'r');

for i = 1:length(signal)
    labels{i} = strcat(signal{i}, ' (', units{i}, ')'); 
end
legend(labels{1,1},labels{1,3});

xlabel('Time (sec)');
grid on;
%% Interpolation to 500Hz
% https://www.google.com.tw/url?sa=t&rct=j&q=&esrc=s&source=web&cd=13&ved=0ahUKEwjnkLL_npfUAhVDoZQKHRsmAjcQFgiAATAM&url=https%3A%2F%2Fmirlab.org%2Fjang%2Fbooks%2FmatlabProgramming4guru%2Fslide%2F09-%25E5%2585%25A7%25E6%258F%2592%25E6%25B3%2595.ppt&usg=AFQjCNF1yMD4McCGM7iTxCQuqiNDQAZo5Q
% http://blog.sina.com.cn/s/blog_4c7482f101009vm2.html

clc;
close all;
clear;
Name ='039m';
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
val3 = interp(val(:,3),4);
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
plot(x(1:range), val(1:range,3),'r');
hold on
plot(x1(1:4*range), val1(1:4*range),'cx');
plot(x1(1:4*range), val3(1:4*range),'y+');
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
%% Preprocessing
figure(1)
range = 750000*4; 
val16000s = val1(1:range);
val36000s = val3(1:range);
x16000s = x1(1:range);
plot(x16000s(1:1000),val16000s(1:1000));
hold on;
plot(x16000s(1:1000),val36000s(1:1000),'r');
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
Bp= val36000s; 
val36000sf = filter(b,1,Bp);
plot(x16000s(1:3000),val16000s(1:3000));
hold on;
plot(x16000s(1:3000),val36000s(1:3000),'r');
hold on;
plot(x16000s(1:3000),val16000sf(1:3000),'c');
hold on;
plot(x16000s(1:3000),val36000sf(1:3000),'y');
grid on;
xlabel('Time (sec)');
ylabel('Amp (mv)');

%% find beats

% Method 1:
Fs=500;
for i=1:30
    so_chan(val16000sf(1+(i-1)*100000:i*100000),Fs) 
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
% Method 2:





% 
% [R, Rt] = findpeaks(val16000sf);                 % Find R-Waves & Times
% % EKGf=(fix(val16000sf));
% % [R, Rt] = findpeaks(EKGf, Fs, 'MinPeakHeight',500);
% figure(3)                                                           % Plot Filtered EKG
% plot(x16000s(1:3000),val16000sf(1:3000),'y');

% grid on
% % display('***This example will write a  Ex1.dat and Ex1.hea file to your current directory!')
% s=input('Hit "ctrl + c" to quit or "Enter" to continue!');
% 
% %Generate 3 different signals and convert them to signed 16 bit in WFDB format
% clear all;clc;close all
% N=1024;
% Fs=48000;
% tm=[0:1/Fs:(N-1)/Fs]';
% adu='V/mV/V';
% info='Example 1';
% 
% 
% %First signal a ramp with 2^16 unique levels and is set to (+-) 2^15 (Volts)
% %Thus the header file should have one quant step equal to (2^15-(-2^15))/(2^16) V.
% sig1=double(int16(linspace(-2^15,2^15,N)'));
% 
% %Second signal is a sine wave with 2^8 unique levels and set to (+-) 1 (mV)
% %Thus the header file should one quant step equal a (1--1)/(2^8)  adu step
% sig2=double(int8(sin(2*pi*tm*1000).*(2^7)))./(2^7);
% 
% %Third signal is a random binary signal set to to (+-) 1 (V) with DC (to be discarded)
% %Thus the header file should have one quant step equal a 1/(2^15) adu step.
% sig3=(rand(N,1) > 0.97)*2 -1 + 2^16;
% 
% %Concatenate all signals and convert to WFDB format with default 16 bits (empty brackets)
% sig=[sig1 sig2 sig3];
% mat2wfdb(sig,'Ex1',Fs,[],adu,info)
% 
% % %NOTE: If you have WFDB installed you can check the conversion by
% % %uncomenting and this section and running (notice that all signals are scaled
% % %to unit amplitude during conversion, with the header files keeping the gain info):
% 
% !rdsamp -r Ex1 > foo
% x=dlmread('foo');
% subplot(211)
% plot(sig)
% subplot(212)
% plot(x(:,1),x(:,2));hold on;plot(x(:,1),x(:,3),'k');plot(x(:,1),x(:,4),'r')

%%
Fs=500;
 i=30
    so_chan(val16000sf(1+(i-1)*100000:i*100000),Fs) 
    ansy{i} = ans{:};

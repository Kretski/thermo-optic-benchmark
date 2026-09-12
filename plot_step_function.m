clear;

% file1 is for the optical response; 
filename1='STEP-20DB-V0.4-8.csv';
sdata=csvread(filename1,2,0);
x1=sdata(:,1)*1e6+0.8;  % time axis, to shift to zero 
y1=sdata(:,2); %optical transmission 
y1n=(y1-min(y1))/(max(y1)-min(y1))-0.03; %adjust for the bias

%file2 is for the electrical response;
%the two channels in the osilloscope is not synchronized 
filename2='STEP-20DB-V0.4-IN8.csv';
sdata=csvread(filename2,2,0);
x2=sdata(:,1)*1e6+553.14; % time axis, to shift to zero (not triggered)
y2=sdata(:,2); 
y2n=(y2-min(y2))/(max(y2)-min(y2))-0.07; %adjust for bias

figure; 
plot(x1,y1n/max(y1n)*1.02);
hold on; 
plot(x2,y2n/max(y2n))

xlim([-6 21]); % unit: us
ylim([-0.05 1.01])
xlabel('Time (\mu s)');
ylabel('Normalized response')

% backl of envelope calculation for filter choices for simultaneous 2p and
% opto
%
% Dec 2017 jvoigts@mit.edu


% Three spectra to consider:
% GCamP
% mRuby/alexa
% JAWs


%% make Gcamp and Rcamp spectra
I= imread('gcamp_spectrum.jpg');
figure(1); clf; 
plot(350,0);
hold on;
image([282,1202],[-0.26,1.08],flipud(I));


wl=[300,350,450,480, 490,500,505,510,512,515,518,520,525,532,538,543,550,560,572,590,620,650,700,900];
a=[0,0,0,0.01,0.1,0.4,0.7,0.95,1,1,0.9,0.8,0.6,0.45,0.38,0.35,0.32,0.2,0.1,0.05,0.01,0.005,0,0];

plot(wl,a,'b-','LineWidth',2);

spectra(1).wl=wl;
spectra(1).a=a;
spectra(1).label='gcamp';


wl=[300,530,550,560,570,580,585,590,594,596,598,601,608,615,627,638,652,670,690,712,730,740,780,800,900];
a=[0,0,0.01,0.06,0.23,0.6,0.8,0.95,1,1,0.99,0.95,0.8,0.65,0.5,0.4,0.3,0.18,0.1,0.05,0.025,0.02,0.005,0,0];
plot(wl,a,'r-','LineWidth',2);

spectra(2).wl=wl;
spectra(2).a=a;
spectra(2).label='rcamp';

%% make JAWs spectrum
I= imread('jaws_spectrum.png');
figure(1); clf; 
plot(350,0);
hold on;
image([373,730],[-2.3,10.7],flipud(I));


wl=[300,350,400,430,460,480,495,515,530,545,560,580,590,605,615,630,645,660,670,680,690,700,900];
a=[0,2,2.7,3,4,5,5.5,6,7,8,9,9.9,10,9.5,8.7,7,4,2,1.1,0.5,0.2,0.1,0];


plot(wl,a,'r-','LineWidth',2);
a=a./max(a);

spectra(3).wl=wl;
spectra(3).a=a;
spectra(3).label='jaws';

%% load filter spectra
filternames={'BLP01-633R','FF01-575_59','FF01-612_SP','FF01-640_20','FF01-640_40','FF611-SDi01'};

for i=1:numel(filternames)
filename = [filternames{i},'_Spectrum.txt'];
delimiter = '\t';
startRow = 5;
formatSpec = '%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

if (min(dataArray{1})>300)
    dataArray{1}=[300;dataArray{1}]';
    dataArray{2}=[dataArray{2}(1);dataArray{2}]';
end;
filters(i).wl=dataArray{1};
filters(i).a=dataArray{2};
filters(i).label=filternames{i}; 
end;

%% do calculations.
figure(2); clf; 
subplot(221); hold on; grid on;
title('JAWs excitation, cleanup+dichroic');

wl=linspace(300,900,100);

a_jaws=interp1(spectra(3).wl,spectra(3).a,wl);
f_cleanup=5;
t_cleanup=interp1(filters(f_cleanup).wl,filters(f_cleanup).a,wl);
plot(wl,a_jaws,'r');

plot(wl,t_cleanup,'b');

f_dichroic=6;
t_dichroic=interp1(filters(f_dichroic).wl,filters(f_dichroic).a,wl);
jaws_e=a_jaws.*t_cleanup.*(1-t_dichroic);
jaws_e(isnan(jaws_e))=0;
plot(wl,t_dichroic,'b--');

text(700,0.5,num2str(sum(jaws_e)./sum(a_jaws)));

f_block=3;
t_block=interp1(filters(f_block).wl,filters(f_block).a,wl).^2;
plot(wl,t_block,'b');
plot(wl,jaws_e,'k','LineWidth',1.5);
legend('JAWs spectrum','cleanup filter','dichroic','blocking filter','JAWs efficiency');

subplot(222); 
semilogy(wl,t_cleanup.*t_dichroic.*t_block,'k');
hold on;
semilogy(wl,t_cleanup,'r');
title('OD clock of light source to PMT');
grid on;
legend('total blocking','cleanup filter');


subplot(223);  hold on;
title('gcamp trough blocking & dichroic');
a_gcamp=interp1(spectra(1).wl,spectra(1).a,wl);
plot(wl,a_gcamp,'r');
gcamp_e=a_gcamp.*t_block.*t_dichroic;
gcamp_e(isnan(gcamp_e))=0;
plot(wl,gcamp_e,'k','LineWidth',1.5);
plot(wl,t_block,'b');
plot(wl,t_dichroic,'b--');
text(700,0.5,num2str(sum(gcamp_e)./sum(a_gcamp)));

legend('GCaMP spectrum','GCaMP efficiency','blocking filter','dichroic');

subplot(224);  hold on;
title('rcamp trough blocking & dichroic');
a_rcamp=interp1(spectra(2).wl,spectra(2).a,wl);
plot(wl,a_rcamp,'r');
rcamp_e=a_rcamp.*t_block.*t_dichroic;
rcamp_e(isnan(rcamp_e))=0;
plot(wl,rcamp_e,'k','LineWidth',1.5);
plot(wl,t_block,'b');
plot(wl,t_dichroic,'b--');
text(700,0.5,num2str(sum(rcamp_e)./sum(a_rcamp)));
plot(wl,jaws_e,'r','LineWidth',1.5);
legend('RCaMP spectrum','RCaMP efficiency','blocking filter','dichroic','JAWs efficiency');

saveas(gcf,'filter_overview.png')






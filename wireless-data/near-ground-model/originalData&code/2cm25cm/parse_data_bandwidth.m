clear all
clc

 
OrgBand=5e9;  %Orignal Bandwidth from 3Ghz to 8Ghz
Band=5e9;    %Slected Bandwidth
Centralband=0e9; %central band
B_end=8e9; %end band

if Centralband~=0 
   B_start=Centralband -Band/2; 
   pathName = sprintf('centralf%dHz',Centralband);
   logf=Centralband;
elseif B_end~=0e9
   B_start=B_end - Band;
   pathName = sprintf('endf%dHz',B_end);
   logf=B_end;
else
   B_start=3e9; %Low frequency of select Bandwidth
   pathName = sprintf('startf%dHz',B_start);
   logf=B_start;
end

noi = -90;   %noise threshold
side =-20;
%measurement distance by dsh
real_distance =6; %unit: m
distance_bias = 0.4;
height = 0.04;       %unit: m
title_time = sprintf('%s distance:%2fm, height: %2fm,start f:%2fHz, \nbandwidth: %2fHz','Time Domain-',real_distance,height,B_start,Band);
fileName = sprintf('%s_bandwidth%dHz',pathName,Band);
c1 = 299792458;


tstart=0e-9;
if Band>0.3e9
    tstop = 30e-9;
else
    tstop = 100e-9;
end

% n_tstart=200e-9;
% n_tstop=300e-9;


 secPeak=1.62*10^(-9);
 peak_width=1;
 flag_fig = 1;
 ampResult = [];
 delayResult = [];
 index = [];
 ftoa=[];
 ftoa_delay=[   ];
 ftoa_amp=[    ];
 
 TOA_dis=[   ];
 
 firstPeakDelay = [    ];
 firstPeakAmp = [   ];
 


 
 P_num=fix((Band/OrgBand)*1601);
 if B_start~=3e9
    P_start=fix(((B_start-3e9)/OrgBand)*1601);
 else
    P_start=1;
 end
 
 P_stop=P_start+P_num-1;
 
 WaveLength= 3e8/Band;
 
 RSS=[];
 PK=[];
 PKgain=[];
 PKdis=[];
 noise_PKgain=[];
 noise_PKgain_1=[];
 
 maxPKgain=[];
 
 DDP=[];
 NDDP=[];
 UDP=[];
 DDP_num=0;
 NDDP_num=0;
 UDP_num=0;

 Fig=1;
 number=500;
 bias=0;
 Profile_number=min(10,number);
 

 
 noise_thred = -90;  %noise threshold
 side =-20;
 
if Fig==1
 h1=figure(1);hold on; grid on;
 xlabel('Delay (s)');
 ylabel('linear Path Loss (w)');
 title(title_time);

 h5=figure(5);hold on;grid on;
 xlabel('Delay (s)');
 ylabel('Path Loss (dB)');
 title(title_time);
 
 h6=figure(6);hold on;grid on;  % 用来保证线性化坐标是一致的。
 xlabel('Delay (s)');
 ylabel('Linear Path Loss (w)');
 title(title_time);
end

log_amp = zeros(1,500);

for i=1:number

 fname = ['scen3_pt' num2str(i) '.s1p'];
 
%  [Hf1, f1] = load_chmeas_s1p_dB( fname, flag_fig );
 
   %%%%%%%%%%%%%Jie He(db)find RSS %%%%%%%%%%%%%%%%%%%%%
   
   fid = fopen(fname, 'rt');
   
   if fid == -1 
    disp(['File cannot be opened !']);
    Hf = 0;  f = 0;
    return;
   end;
   
   while( 1 )
    temp_str = fgetl(fid);  % read in a line of text.
    
    if temp_str(1) == '!'
        if flag_fig == 1
            disp(temp_str);
        end;
    else 
        if temp_str(1) == '#'
            tmp_data = fscanf(fid, '%g %g %g', [3 inf] );
            fclose(fid);
            
            tmp_data = tmp_data.';
            f = tmp_data(:,1);
            amp = 10.^(tmp_data(:,2)/20);
      % Channel Transfer Function measured by VNA
            %      Hf = 10.^(tmp_data(:,2)/20).*exp(1j*tmp_data(:,3)*pi/180);
            Hf = amp.*exp(1j*tmp_data(:,3)*pi/180);
            break;    
        else 
            if feof(fid)
                fclose(fid);
                Hf = 0;  f = 0;
                return;
            end;
        end;
    end;
end;
  
 f_dB=20*log10(abs(amp))-bias;
 RSS_dB=mean(f_dB);
 
 RSS=[RSS,RSS_dB];
  
  
  %%%%%%%%%%%%%Jie He(db)find RSS %%%%%%%%%%%%%%%%%%%%%
  

 %[zt_han, t] = czt_hanning( f1, Hf1, tstart, tstop, 1, 1601*1000);

%  [zt_han, t] = czt_hanning( f1, Hf1, tstart, tstop, 1, 1601);

 Hf=Hf(P_start:P_stop);
 f=f(P_start:P_stop);
 [zt_han, t] = czt_hanning( f, Hf, tstart, tstop, 1, 1601*10);
 
 time_dB = 20*log10(abs(zt_han))-bias;

 %%%%%%%%%%%%%Jie He(db)find noise %%%%%%%%%%%%%%%%%%%%%

%  Num_start=fix((n_tstart/tstop)*length(t));
%  Num_stop=fix((n_tstop/tstop)*length(t));
%  
% 
%  
%  noise_index = pkd_cir(time_dB( Num_start:Num_stop), noise_thred, side, peak_width)+Num_start;
%  
%  if length(noise_index)~=0
%      noise_PKgain_1=[];
%      for k=1:length(noise_index)
%         noise_PKgain_1=[noise_PKgain_1,20*log10(abs(zt_han(noise_index(k))))-bias];
%         noise_PKgain=[noise_PKgain, 20*log10(abs(zt_han(noise_index(k))))-bias];
%      end
%  end
 %%%%%%%%%%%%%Jie He(db)find noise %%%%%%%%%%%%%%%%%%%%%


 %%%%%%%%%%%%%Jie He(db) finde peak %%%%%%%%%%%%%%%%%%%%%

 
 index = pkd_cir(time_dB, noi, side, peak_width);
 

 %%%%%%%%%%%%%orignal(mw)%%%%%%%%%%%%%%%%%%%%% 
%  index = pkd_cir(abs(zt_han), noi, side, peak_width);
 if index == 0 
     continue
 end
        

 
 ftoa_delay = [ftoa_delay t(index(1))];
 ftoa_amp = [ftoa_amp 20*log10(abs(zt_han(index(1))))-bias];
 
 % Plot Time Response in Time Domain - disable for looping
 direct_delay = real_distance/ c1
 bias_delay = ftoa_delay(1) - direct_delay
 if Fig==1 && Profile_number~=0
     figure(5);hold on;grid on;
     plot(t,10.^(time_dB./10),'b');
     figure(6);hold on;grid on;
     plot(t,10.^(time_dB./10),'b');
     figure(1);hold on;grid on;
     plot(t,time_dB,'b');
     
     figure(5);hold on;
     plot(t(index(1:length(index))),10.^((20*log10(abs(zt_han(index(1:length(index)))))-bias)./10),'bo');
     plot(ftoa_delay,10.^(ftoa_amp./10),'ro');
     
     figure(6);hold on;
     plot(t(index(1:length(index))),10.^((20*log10(abs(zt_han(index(1:length(index)))))-bias)./10),'bo');
     plot(ftoa_delay,10.^(ftoa_amp./10),'ro');
     
     figure(1);hold on;
     plot(t(index(1:length(index))),(20*log10(abs(zt_han(index(1:length(index)))))-bias),'bo');
     plot(ftoa_delay,ftoa_amp,'ro');
     
     t2=t(index(1:length(index)));
     t1=20*log10(abs(zt_han(index(1:length(index)))))-bias;
     str1=sprintf('other paths:\n');
     if length(t2)>=2
         for index1=2: length(t2)
              str1=sprintf('%s time:%e, distance:%f\n',str1,t2(index1)+bias_delay, (t2(index1)+bias_delay)*c1);
              ftoa_delay_ground1 = t2(index1);
              ftoa_amp_ground1 = t1(index1);
              figure(5);hold on;
              plot(ftoa_delay_ground1,10.^(ftoa_amp_ground1./10),'go');
              figure(1);hold on;
              plot(ftoa_delay_ground1,ftoa_amp_ground1,'go');
         end
         secTime = 2*sqrt((real_distance/2 )^2 +height^2)/c1;
         min=t2(2)-secTime;
         ftoa_delay_ground = t2(2);
         ftoa_amp_ground = t1(2);
         for index1 =2:length(t2)
             if abs(t2(index1)-secTime)<min
                  ftoa_delay_ground = t2(index1);
                  ftoa_amp_ground = t1(index1);
             end
         end
         figure(5);hold on;
         plot(ftoa_delay_ground,10.^(ftoa_amp_ground./10),'go');
         figure(6);hold on;
         plot(ftoa_delay_ground,10.^(ftoa_amp_ground./10),'go');
         figure(1);hold on;
         plot(ftoa_delay_ground,ftoa_amp_ground,'go');
     end
     
%      
%      if length(noise_index)~=0
%          plot(t(noise_index),noise_PKgain_1,'bo');
%      end
     
     Profile_number=Profile_number-1;
 end
 
 for k=1:length(index)
    PKgain(i,k)= 20*log10(abs(zt_han(index(k))))-bias;
    PKdis(i,k)=t(index(k))*2.99792458*10^8;
     
    PK(i,2*(k-1)+2) = 20*log10(abs(zt_han(index(k))))-bias;
    PK(i,2*(k-1)+1)=t(index(k))*2.99792458*10^8;
 end
 
 maxPKgain(i)=max(PKgain(i,1:k));
 
 if PKdis(i,1)>=(real_distance + distance_bias)-3e8/Band && PKdis(i,1)<=(real_distance + distance_bias)+3e8/Band 
     if PKgain(i,1)==maxPKgain(i);
         DDP=[DDP i];
         DDP_num=DDP_num+1;
     else  
         NDDP=[NDDP i];
         NDDP_num=NDDP_num+1;
     end  
 else
     UDP=[UDP i];
     UDP_num=UDP_num+1;
 end
 

 
 
 
%  for k=1:length(index)
%  ftoa(k)= min(t(index(k)));
%  gain = abs(zt_han(index(k)));
%  power=10*log10(sum(gain.*gain));s
%  end
%  
%  for ii=1:length(index)
%      pk_gain = abs(zt_han(index(ii)));
%      pk_delay = t(index(ii));
%      fprintf(Rfid,'%g                %g\n\n',pk_gain,pk_delay);
%  end
 
end



ftoa_dist=ftoa_delay*2.99792458*10^8;
TOA_dis=sort(ftoa_dist);
TOA_Error=ftoa_dist-(real_distance + distance_bias);
locate=find(TOA_Error<0);
TOA_Error(locate)=[];
% locate=find(TOA_Error>0.5);
% TOA_Error(locate)=[];

Mean_RSS=mean(RSS)
Mean_FP=mean(ftoa_amp)
Mean_dis=mean(TOA_Error)
Mean_error=mean(TOA_Error);

Var_RSS=var(RSS)
Var_FP=var(ftoa_amp)
Var_dis=var(TOA_dis)
Var_error=var(TOA_Error)

% log the file of amplitude
pathFile= sprintf('..\\%s',pathName);
if exist(pathFile)==0   
    pathFile1=sprintf('mkdir %s',pathFile);
    system(pathFile1);
end
pathFile= sprintf('..\\%s\\amptitude.txt',pathName);
fid = fopen(pathFile, 'a+');

fprintf(fid, '%f %f %d %d ', height,real_distance,logf,Band );

fprintf(fid, '%f ', ftoa_amp);
fprintf(fid,'%f %f', Mean_FP,Var_FP);
fprintf(fid, '\n');
fclose(fid);

% log the file of ranging error
% log the file of amplitude
pathFile= sprintf('..\\%s',pathName);
if exist(pathFile)==0   
    pathFile1=sprintf('mkdir %s',pathFile);
    system(pathFile1);
end
pathFile= sprintf('..\\%s\\ranging_error.txt',pathName);
fid1= fopen(pathFile, 'a+');

fprintf(fid1, '%f %f %d %d ', height,real_distance,logf,Band );
fprintf(fid1, '%f ', TOA_Error);
fprintf(fid,'%f %f', Mean_error,Var_error);
fprintf(fid1, '\n');
fclose(fid1);

% log the file of ranging error
% log the file of amplitude
pathFile= sprintf('..\\%s',pathName);
if exist(pathFile)==0   
    pathFile1=sprintf('mkdir %s',pathFile);
    system(pathFile1);
end
pathFile= sprintf('..\\%s\\average_summary.txt',pathName);
fid2 = fopen(pathFile, 'a+');
fprintf(fid2, '%f %f %d %d ', height,real_distance,logf,Band );
fprintf(fid2, '%f %f %f %f', Mean_FP,Var_FP,Mean_error,Var_error);
fprintf(fid2, '\n');
fclose(fid2);

% if length(noise_PKgain)~=0
%     Mean_Noise=mean(noise_PKgain)
% end
% Max_Noise=max(noise_PKgain);




%  figure(4);hold on;grid on;
%  title(' First Path Path Loss versus TOA distance');
%  xlabel('TOA distance (m)');
%  ylabel('Path Loss (dB)');
% 
%  ftoa_dist=ftoa_delay*2.99792458*10^8;
%  figure(4);hold on;
%  plot(ftoa_dist,ftoa_amp,'*');

 
 
 
 m=1:1:length(TOA_Error);
 
 h3=figure(3); hold on; grid on;
 title('TOA distance in sequence');
 xlabel('Sequence');
 ylabel('TOA distance(m)');
 plot(m,TOA_Error,'*-');
 
 strange=[];
 sm=0;
 for i=1:length(ftoa_dist)
     
     if ftoa_dist(i)<(real_distance + distance_bias)
         sm=sm+1;
         strange(1,sm)= i;
         strange(2,sm)=ftoa_dist(i);
     end
 end
  
 
%  if(sm > 0)
%      figure(4);hold on;grid on;
%      xlabel('Delay (s)');
%      ylabel('Path Loss (dB)');
%      title('Time Domain');
%      for i=1:sm
%          fname = ['scen3_pt' num2str(strange(1,i)) '_1.s1p'];
%          [Hf1, f1] = load_chmeas_s1p_dB( fname, flag_fig );
%          [zt_han, t] = czt_hanning( f1, Hf1, tstart, tstop, 1, 1601);
%          time_dB = 20*log10(abs(zt_han))-bias;
%          figure(4);hold on;grid on;
%          plot(t,time_dB);
% 
%           %%%%%%%%%%%%%Jie He(db)%%%%%%%%%%%%%%%%%%%%%
% 
%          index = pkd_cir(time_dB, noi, side, peak_width);
% 
%          if index == 0 
%              continue
%          end
% 
% 
%          figure(4);hold on;
%          plot(t(index(1:length(index))),20*log10(abs(zt_han(index(1:length(index)))))-bias,'bo');
%          plot( t(index(1)),20*log10(abs(zt_han(index(1))))-bias,'ro'); 
%      end
%  end
 






%   for k=1:length(index)
%   ftoa(k)= min(t(index(k)));
%   gain = abs(zt_han(index(k)));
%   power=10*log10(sum(gain.*gain));
%   end
  

%    for ii=1:length(index)
%       pk_gain = abs(zt_han(index(ii)));
%       pk_delay = t(index(ii));
%   fprintf(Rfid,'%e                %e\n\n',pk_gain,pk_delay);
%   end
%   fclose(Rfid);  
 
  
 % Separate Tap Amplitude and Delay
%  Sfid = fopen('Result.txt','r');
%  tResult = fscanf(Sfid,'%g');
%  for tNum = 1:length(tResult)
%      test = tResult(tNum);
%      if mod(tNum,2)==0
%      delayResult = [delayResult test];          % Save Delay info to Matrix delayResult[]
%      else ampResult = [ampResult test];         % Save Amplitude info to Matrix ampResult[]
%      end
%  end
 
%  fclose(Rfid);
 

%  for tPeak=1:length(delayResult)
%      if(delayResult(tPeak)>=0 && delayResult(tPeak)<=secPeak)
%          firstPeakDelay = [firstPeakDelay delayResult(tPeak)];
%          firstPeakAmp = [firstPeakAmp ampResult(tPeak)];
%      end
%  end
 
%  figure(1); hold on;grid on;
%  Amp_dB = 20*log10(abs(ampResult));
% % plot(delayResult,Amp_dB,'*');
%  title('Path Amplitude variance versus Path Delay With Water in Phantom');
%  xlabel('Path Delay (ns)');
%  ylabel('Path Ampiltude (dB)');
%  %axis([1.1e-9 1.6e-9 -80 -35]);
 
%  actDist=0.2286; % thickness of phantom 0.2286m
%  figure(3);hold on;grid on;
%  firstPeakAmp_dB = 20*log10(abs(firstPeakAmp));
%  firstPeakDist = firstPeakDelay*3*(10^9)-actDist;
%  plot(firstPeakDist,firstPeakAmp_dB,'*');
%  title(' Distance Measurement Error versus Path Delay With Water in Phantom');
%  xlabel('Path Delay (m)');
%  ylabel('Path Ampiltude (dB)');
% axis([1.35e-9 1.4e-9 -80 -65]);


%  ftoa_dist=ftoa_delay*3*10^8;
figure(5)
ymax=get(gca,'ylim');
clear min
direct_delay = real_distance/ c1
refract_delay = ((real_distance/2)^2+height^2)^0.5 * 2 /c1
bias_delay = ftoa_delay(1) - direct_delay
direct_delay = direct_delay + bias_delay
refract_delay = refract_delay+ bias_delay
scope = ymax(1):(ymax(2)-ymax(1))/100:ymax(2)
delay = ones(size(scope))* direct_delay(1)
plot(delay, scope,'r')
delay = ones(size(scope))* refract_delay(1)
plot(delay, scope,'b')
xlim([1.5e-8,3.0e-8])
str = sprintf('direct delay - refection delay :%e\naverage ranging error: %f', refract_delay -direct_delay,Mean_dis);
str = sprintf('%s \n %s',str,str1);
annotation('textbox', [.5 .8, .1, .1], 'string', str, 'fontsize',9, ...
     'color', 'm', 'edgecolor', 'none');
 
figure(6)
ymax=get(gca,'ylim');
clear min
direct_delay = real_distance/ c1
refract_delay = ((real_distance/2)^2+height^2)^0.5 * 2 /c1
bias_delay = ftoa_delay(1) - direct_delay
direct_delay = direct_delay + bias_delay
refract_delay = refract_delay+ bias_delay
scope = ymax(1):(ymax(2)-ymax(1))/100:ymax(2)
delay = ones(size(scope))* direct_delay(1)
plot(delay, scope,'r')
delay = ones(size(scope))* refract_delay(1)
plot(delay, scope,'b')


figure(1)
ymax=get(gca,'ylim');
clear min
direct_delay = real_distance/ c1
refract_delay = ((real_distance/2)^2+height^2)^0.5 * 2 /c1
bias_delay = ftoa_delay(1) - direct_delay
direct_delay = direct_delay + bias_delay
refract_delay = refract_delay+ bias_delay
scope = ymax(1):(ymax(2)-ymax(1))/100:ymax(2)
delay = ones(size(scope))* direct_delay(1)
plot(delay, scope,'r');
delay = ones(size(scope))* refract_delay(1);
plot(delay, scope,'b');
% save fig file
pathFile= sprintf('.\\%s',pathName);
if exist(pathFile)==0   
    pathFile1=sprintf('mkdir %s',pathFile);
    system(pathFile1);
end


fileNameofFig1=sprintf('%s\\%s_h%d_dis%dfig1.fig',pathFile,fileName,height*100,real_distance*100);

saveas(h1,fileNameofFig1);
 fileNameofFig3=sprintf('%s\\%s_h%d_dis%dfig3.fig',pathFile,fileName,height*100,real_distance*100);


saveas(h3,fileNameofFig3);
fileNameofFig5=sprintf('%s\\%s_h%d_dis%dfig5',pathFile,fileName,height*100,real_distance*100);
saveas(h5,fileNameofFig5);
fileNameofFig6=sprintf('%s\\%s_h%d_dis%dfig6',pathFile,fileName,height*100,real_distance*100);
saveas(h6,fileNameofFig6);


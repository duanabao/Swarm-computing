clear all
clc

 
OrgBand=5e9;  %Orignal Bandwidth from 3Ghz to 8Ghz
Band=5e9;    %Slected Bandwidth
Centralband=5e9; %central band
B_end=0e9; %end band

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
side =-32;
%measurement distance by dsh
real_distance =0.5; %unit: m
distance_bias = 0.4;
height = 0.02;       %unit: m
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
 path=[];
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
 delayofFirsttwopath = zeros(1,number);

 
 noise_thred = -90;  %noise threshold
 side =-32;
 
 
if Fig==1
 h1=figure(1);hold on; grid on;
 xlabel('Delay (s)');
 ylabel('linear Path Loss (dB)');
 title(title_time);

 h5=figure(5);hold on;grid on;
 xlabel('Delay (s)');
 ylabel('Path Loss (w)');
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

 %[zt_han, t] = czt_hanning( f, Hf, tstart, tstop, 1, 1601);

% Hf=Hf(P_start:P_stop);
% f=f(P_start:P_stop);
 [zt_han, t] = czt_hanning( f, Hf, tstart, tstop, 1, 1601*10);
 
 time_dB = 20*log10(abs(zt_han))-bias;

 %%%%%%%%%%%%%Jie He(db)find noise %%%%%%%%%%%%%%%%%%%%%

 % Num_start=fix((n_tstart/tstop)*length(t));
 % Num_stop=fix((n_tstop/tstop)*length(t));
  
 
  
 % noise_index = pkd_cir(time_dB( Num_start:Num_stop), noise_thred, side, peak_width)+Num_start;
  
 % if length(noise_index)~=0
  %    noise_PKgain_1=[];
   %   for k=1:length(noise_index)
    %     noise_PKgain_1=[noise_PKgain_1,20*log10(abs(zt_han(noise_index(k))))-bias];
     %    noise_PKgain=[noise_PKgain, 20*log10(abs(zt_han(noise_index(k))))-bias];
      %end
 %end
 %%%%%%%%%%%%%Jie He(db)find noise %%%%%%%%%%%%%%%%%%%%%


 %%%%%%%%%%%%%Jie He(db) finde peak %%%%%%%%%%%%%%%%%%%%%

 
index = pkd_cir(time_dB, noi, side, peak_width);
 

 %%%%%%%%%%%%%orignal(mw)%%%%%%%%%%%%%%%%%%%%% 
 % index = pkd_cir(abs(zt_han), noise, side, peak_width);
 if index == 0 
     continue
 end
        

 for k=1:length(index)
    PKgain(i,k)= 20*log10(abs(zt_han(index(k))))-bias;
 for j=1:k
[ftoa_amp(i),j]=max(PKgain(i,1:k));

 end
 end

     
 delay_amp(i,1:k-j+1)=PKgain(i,j:k);
 path=[path k-j+1];
 delay(i,1:k-j+1)=t(index(j:k));
  
 
 ftoa_delay = [ftoa_delay t(index(j))];

 
 
 
 % Plot Time Response in Time Domain - disable for looping
 direct_delay = real_distance/ c1
 bias_delay = ftoa_delay(1) - direct_delay
 t2=t(index(j:length(index)));
 t1=20*log10(abs(zt_han(index(j:length(index)))))-bias;
 if length(t2)>=2
         delayofFirsttwopath(i)=t2(2)-t2(1);
        
end
 if Fig==1 && Profile_number~=0
     figure(5);hold on;grid on;
     plot(t,10.^(time_dB./10),'b');
     figure(6);hold on;grid on;
     plot(t,10.^(time_dB./10),'b');
     figure(1);hold on;grid on;
     plot(t,time_dB,'b');
     
     figure(5);hold on;
     plot(t(index(j:length(index))),10.^((20*log10(abs(zt_han(index(j:length(index)))))-bias)./10),'bo');
     plot(ftoa_delay,10.^(ftoa_amp./10),'ro');
     
     figure(6);hold on;
     plot(t(index(j:length(index))),10.^((20*log10(abs(zt_han(index(j:length(index)))))-bias)./10),'bo');
     plot(ftoa_delay,10.^(ftoa_amp./10),'ro');
     
     figure(1);hold on;
     plot(t(index(j:length(index))),(20*log10(abs(zt_han(index(j:length(index)))))-bias),'bo');
     plot(ftoa_delay,ftoa_amp,'ro');
     
     t2=t(index(j:length(index)));
     t1=20*log10(abs(zt_han(index(j:length(index)))))-bias;
     str1=sprintf('other paths:\n');
     if length(t2)>=2
         delayofFirsttwopath(i)=t2(2)-t2(1);
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
  %    if length(noise_index)~=0
   %       plot(t(noise_index),noise_PKgain_1,'bo');
    % end
     
     Profile_number=Profile_number-1;
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
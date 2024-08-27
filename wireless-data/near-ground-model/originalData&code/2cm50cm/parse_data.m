clear all
clc

 
OrgBand=5e9;  %Orignal Bandwidth from 3Ghz to 8Ghz
 
B_start=3e9; %Low frequency of select Bandwidth
 
Band=0.1e9;    %Slected Bandwidth

noi = -70;   %noise threshold
side =-20;
%measurement distance by dsh
real_distance = 6;
distance_bias=0.4;
hight = 0.02;
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
 figure(5);hold on;grid on;
 xlabel('Delay (s)');
 ylabel('Path Loss (dB)');
 title('Time Domain');
end



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

 if Fig==1 && Profile_number~=0
     figure(5);hold on;grid on;
     plot(t,time_dB,'b');
     figure(5);hold on;
     
     
     plot(t(index(1:length(index))),20*log10(abs(zt_han(index(1:length(index)))))-bias,'bo');
     plot(ftoa_delay,ftoa_amp,'ro');
     
     t2=t(index(1:length(index)));
     if length(t2)>=2
         t1=20*log10(abs(zt_han(index(1:length(index)))))-bias;
         secTime = 2*sqrt((real_distance/2 )^2 +hight^2)/c1;
         min=t2(2)-secTime;
         ftoa_delay_ground = t2(2);
         ftoa_amp_ground = t1(2);
         for index1 =2:length(t2)
             if abs(t2(index1)-secTime)<min
                  ftoa_delay_ground = t2(index1);
                  ftoa_amp_ground = t1(index1);
             end
         end
     
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
% locate=find(TOA_Error<0);
% TOA_Error(locate)=[];
Mean_RSS=mean(RSS)
Mean_FP=mean(ftoa_amp)
Mean_dis=mean(TOA_dis)
Mean_error=mean(TOA_Error)

Var_RSS=var(RSS)
Var_FP=var(ftoa_amp)
Var_dis=var(TOA_dis)
Var_error=var(TOA_Error)



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
 
 figure(3); hold on; grid on;
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

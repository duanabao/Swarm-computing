%
% This program is used to load 8753D Network Analyzer Measurement data.
% Read S21 data from S1P file.  LogMag/Angle.
%

function [ Hf, f ] = load_chmeas_s1p_dB( fname, flag_fig) 

%fid = fopen(fname, 'rt');
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
    
% %Plot figure in frequency domain - disable while looping
% if flag_fig == 1
%     tmp_f = f*1e-9;
%     mag_dB = 10*log10(abs(Hf));
%     phs = angle(Hf.');
% 
%     figure; hold on; box on;
%    subplot(2,1,1);  plot(tmp_f, mag_dB);
%     subplot(2,1,1);  plot(tmp_f, mag_dB);
%     xlabel('frequency (GHz)');
%     ylabel('Magnitude (dB)');
%     title(fname);
% 
%     subplot(2,1,2);  plot(tmp_f, phs);  
%     xlabel('frequency (GHz)');  
%     ylabel('angle (radian)');
% end;

return;
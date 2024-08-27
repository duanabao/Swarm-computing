%
% Peak detection on channel impulse response.
%
% input: 
%       ht: channel impulse response
%       noi: threshold for noise std
%       side: sidelobe amplitude for window functions
%           Rec: -13dB, Hanning: -32dB, Hamming: -43dB  
%       peak_width: time resolution of peak in units of dt

function [ peak_index ] = pkd_cir(ht, noi, side, peak_width)
% peak_width is not used in this version.

len_t = length(ht);
peak = max(ht);

peak_index = 0;
count = 0;
i = 2;
while(1)
    %%%%%%%%%%%%%Orignal(mw)%%%%%%%%%%%%%%%%%%%%%
%     if ht(i)>ht(i-1) & ht(i)>ht(i+1) & ht(i)>noi & ht(i)/peak > side  
    %%%%%%%%%%%%%Orignal%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%Jie He(db)%%%%%%%%%%%%%%%%%%%%%
    if ht(i)>ht(i-1) & ht(i)>ht(i+1) & ht(i)>noi & ht(i)-peak > side
         %%%%%%%%%%%%%Jie He%%%%%%%%%%%%%%%%%%%%%
         

        if count == 0
            peak_index = i;
            count = 1;
        else
            peak_index = [peak_index, i];
        end;
        i = i + 1;
    else
        i = i + 1;
    end;

    if i > len_t - 1
        break;
    end;
end;

return;

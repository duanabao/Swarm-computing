% Compute Channel Impulse Response from frequency measurement 
% data using Chirp-Z transform with hanning window.
% modified 03/27/02.

function [ zt_han , t ] = czt_hanning(freq, Zf, tstart, tstop, flag, Nt )

%Nf = length(freq);
Nf = length(freq);
df = (freq(Nf)-freq(1))/(Nf-1);

T = 1/df;

if nargin < 6
 % Nt = 1601;
  Nt = 1601;
end;

if flag == 1
  han = hanning(Nf);
%  han = hann(Nt);
  Zf = (45/23)*Zf(:).*han(:);  % 45/23 is to make the Hanning-window time response peak at 1.
%   Zf = Zf(:).*han(:);
end;

dt = (tstop-tstart)/(Nt-1);
w = exp(1j*2*pi*dt/T);
a = exp(1j*2*pi*tstart/T);

zt_han = (1/Nf)*czt(Zf(:), Nt, w, a);

t = linspace(tstart, tstop, Nt);

return;

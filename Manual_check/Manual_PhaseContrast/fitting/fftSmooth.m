function [s]=fftSmooth(signal,fsmoothValue)
%function frdescp from Oufti

[np, nc] = size(signal);
if nc ~=2, error('S must be of size np-by-2.'); end
if np/2 ~= round(np/2)
   signal(end+1,:) = signal(end, :);
   np = np + 1;
end
x = 0:(np-1);
m = ((-1).^x)';
signal(:,1) = m .* signal(:,1);
signal(:,2) = m .* signal(:,2);
signal = signal(:,1) + sqrt(-1)*signal(:,2);
inverseSignal = fft(signal);

%function ifdescp from Oufti

signalPoints = length(inverseSignal);
x = 0:(signalPoints-1);
m = ((-1).^x)';
d = round((signalPoints - fsmoothValue)/2);
inverseSignal(1:d) = 0;
inverseSignal(signalPoints-d+1:signalPoints) = 0;
zz = ifft(inverseSignal);
s(:,1) = real(zz);
s(:,2) = imag(zz);
s(:,1) = m.*s(:,1);
s(:,2) = m.*s(:,2);
end
function pulsed_pink(fname,FS,stim_dur,pulse_dur)

% 'fname' is the outputfilename
% 'FS' is the sample rate
% 'stim_dur' in s is the total duration of the resulting noise
% 'pulse_dur' in s is the duration of a high pulse, no pulse if variable is absent


nsamples = round(FS * stim_dur); % number of samples
y = pinknoise(nsamples); % genereate pink_noise

% high pass (4 passes)
fpass = 20;
bhi = fir1(512,fpass/(FS/2),'high');
y = filtfilt(bhi,1,y);
y = filtfilt(bhi,1,y);

if exist('pulse_dur','var') % pulse
    freq = 1 / (pulse_dur * 2);
    t = linspace(0,stim_dur,nsamples);
    carre = square(2 * pi * freq * t);
    carre = (carre / 2) + 0.5;
    y = y .* carre;
end



y = y ./ max(abs(y));

% figure
% hold on
% plot(t,y,'k')
% plot(t,carre,'r')

audiowrite(fname,y,FS)

end

% SUBFUNCTIONS
function y = pinknoise(N)

% function: y = pinknoise(N) 
% N - number of samples to be returned in row vector
% y - row vector of pink (flicker) noise samples

% The function generates a sequence of pink (flicker) noise samples. 
% Pink noise has equal energy in all octaves (or similar log bundles) of frequency.
% In terms of power at a constant bandwidth, pink noise falls off at 3 dB per octave. 

% difine the length of the vector
% ensure that the M is even
if rem(N,2)
    M = N+1;
else
    M = N;
end

% generate white noise
x = randn(1, M);

% FFT
X = fft(x);

% prepare a vector for 1/f multiplication
NumUniquePts = M/2 + 1;
n = 1:NumUniquePts;
n = sqrt(n);

% multiplicate the left half of the spectrum so the power spectral density
% is proportional to the frequency by factor 1/f, i.e. the
% amplitudes are proportional to 1/sqrt(f)
X(1:NumUniquePts) = X(1:NumUniquePts)./n;

% prepare a right half of the spectrum - a copy of the left one,
% except the DC component and Nyquist frequency - they are unique
X(NumUniquePts+1:M) = real(X(M/2:-1:2)) -1i*imag(X(M/2:-1:2));

% IFFT
y = ifft(X);

% prepare output vector y
y = real(y(1, 1:N));

% ensure unity standard deviation and zero mean value
y = y - mean(y);
yrms = sqrt(mean(y.^2));
y = y/yrms;

end
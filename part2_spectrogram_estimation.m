% The code for getting the spectrograms
clear all;
close all;

Fs = 5000;  % Sampling Frequency
A = 5;  % amplitude of the signal
c = 340;    % speed of the signal (m/s)
v = 50+randi(200)   % moving source speed (m/s)
n = randi(5);
k = randi(10);
t0 = (n+k*0.1)/Fs  % time delay (s)
x = v*t0  % the distance between source and the sensor (m)
f0 = 350+randi(100)    % signal frequency
delta = randi(2); % duration of the signal (s)
nos = delta*Fs + 1; % Number of samples
m = randi(100);

InputType = 1;   % input type 1=sinus, 0=linear chirp


% Emitted Signal
t = 0:1/Fs:delta;
if InputType == 1   % Sinusoidal Input %
    Se = A*sin(2*pi*f0*t);  % emitted signal
elseif InputType == 0  % Linear Chirp Input %
    Se = A*cos(2*pi*f0*t + 2*pi*400*(t.*t)); % emitted signal
end


e = wgn(1,nos,0); % white noise

% Power of the input signal
Pi = sum(Se.^2)/nos;
% Power of noise signal
Pn = sum(e.^2)/nos;

% Signal to Noise Ratio
SNR = Pi/Pn;

% Received Signal
t = (c/(c-v))*(t-t0); 
if InputType == 1   % Sinusoidal Input %
    Sr = A*sin(2*pi*f0*t) + e;
elseif InputType == 0  % Linear Chirp Input %
    Sr = A*cos(2*pi*f0*t + 2*pi*400*(t.*t)) + e;
end


x=Sr;   % to not change our part1 code
M=220;  % window overlap
N=250;  % window size

[X, w] = freqz(x);

% Normalize the signal
x = x.'/max(abs(x));

% Make the signal a row vector
x = x(:).';

% Number of segments (frames) the signal is divided to.
K = floor((nos-M)/(N-M)); 

%STFT
X = zeros(N,K);
S = zeros(N,K);

win = gausswin(N).';    % Gaussian Window

for k=1:K
X(:,k) = x((k-1)*(N-M)+1:k*N - (k-1)*M).*win;
S(:,k) = fft(X(:,k));
end
S = S(1:N/2,:); % STFT is saved in matrix S (N/2)xK

t =(t0:K-1+t0)*nos/(K*10000);

%Frequency points in Hz
f = (0:N/2-1)*Fs/N;  

% Frequency estimation of the signal from STFT
if InputType == 1
    for i=1:(N/2)
        if abs(S(i,10)) == max(abs(S(:, 10)))
            fsin = i*Fs/N;
        end
    end
    
elseif InputType == 0
     for i=1:(N/2)
        if abs(S(i,2)) == max(abs(S(:, 2)))
            fsin = i*Fs/N;
        end
        if abs(S(i,K-1))== max(abs(S(:, K-1)))
            fchirp2 = i*Fs/N;
        end
     end
end


% Speed calculation from the frequency estimation
Vest = c*(fsin - f0)/fsin  % m/s

% Range calculation from the time of arrival and speed
Xest = Vest*t0  % m


%Plot the spectogram
h = figure('Name','STFT - Spectogram');
colormap('jet');

[T,F] = meshgrid(t,f/1000); % f in KHz.
% For spectogram, square magnitude of STFT is used
surface(T,F,10*log10(abs(S.^2) + eps),'EdgeColor','none');
axis tight;
grid on;
title(['Fs: ',num2str(Fs),', Window Length: ', num2str(N),', Overlap: ', num2str(M)]);
xlabel('Time (sec)');
ylabel('Frequency (KHz)');
colorbar('Limits',[-80, 40]);
cbar_handle = findobj(h,'tag','Colorbar');
set(get(cbar_handle,'YLabel'),'String','(dB)','Rotation',0);
zlim([-80 40]);
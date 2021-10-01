% The Code for determining the parameter ranges
close all;
clear all;

InputType = 1;  % input type 1=sinus, 0=linear chirp
Fs = 5000;  % sampling frequency

% default values
delta = 2; % duration of the signal (s)
f0 = 400;   % signal frequency
n = 7;
k = 3;
t0 = (n+k*0.1)/Fs;  % time delay (s)
v = 120;   % moving source speed (m/s)
x = v*t0;  % the distance between source and the sensor (m)
m = 0;  % Frequency parameter for linear chirp

NT = 10000; % Number of trials
MC = zeros(2, NT);
SNRM = zeros(1, NT);

for i = 1:NT

    % Speed range
    v = 50+randi(200);   % moving source speed (m/s)

    % time delay range
    n = randi(5);
    k = randi(10);
    t0 = (n+k*0.1)/Fs;  % time delay (s)
    
    x = v*t0;  % the distance between source and the sensor (m)
    
    % f0 range
    f0 = 350+randi(100);

    % delta range
    delta = randi(2); % duration of the signal (s)
    
    % m range
    m = randi(100);
    
    [Vest, Xest, SNR] = myestimate(InputType, f0, delta, t0, v, m);
    
    
    MC(1, i) = abs(Vest-v); % Speed estimation errors
    MC(2, i) = abs(Xest-x); % Range estimation errors

%     SNRM(1, i) = SNR;
end

figure
histogram(MC(1,:));
figure
histogram(MC(2,:));
% figure 
% histogram(SNRM);

function [Vest, Xest, SNR] = myestimate(InputType, f0, delta, t0, v, m)
    
    Fs = 5000;  % Sampling Frequency
    nos = delta*Fs + 1; % Number of samples
    A = 5;  % amplitude of the signal
    c = 340;    % speed of the signal (m/s)
    
    % Emitted signal
    t = 0:1/Fs:delta;
    if InputType == 1   % Sinusoidal Input %
     Se = A*sin(2*pi*f0*t);  % emitted signal
    elseif InputType == 0  % Linear Chirp Input %
        Se = A*cos(2*pi*f0*t + 2*pi*(t.*t)*m/(2*delta)); % emitted signal
    end
    
    e = wgn(1,nos,0); % white noise
    
    % Power of the input signal
    Pi = sum(Se.^2)/nos;
    % Power of noise signal
    Pn = sum(e.^2)/nos;

    % Signal to Noise Ratio
    SNR = Pi/Pn;
    
    
    % Received signal
    t = (c/(c-v))*(t-t0); 
    if InputType == 1   % Sinusoidal Input %
        Sr = A*sin(2*pi*f0*t) + e;
    elseif InputType == 0  % Linear Chirp Input %
        Sr = A*cos(2*pi*f0*t + 2*pi*(t.*t)*m/(2*delta)) + e;
    end
    
    % STFT Calculation
    
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
            if abs(S(i,1)) == max(abs(S(:, 1)))
                fsin = i*Fs/N;
            end
            if abs(S(i,K-1))== max(abs(S(:, K-1)))
                fchirp2 = i*Fs/N;
            end
         end
    end

    % Speed calculation from the frequency estimation
    Vest = c*(fsin - f0)/fsin;  % m/s

    % Range calculation from the time of arrival and speed
    Xest = Vest*t0;  % m

end

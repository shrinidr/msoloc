function [] = input_simul()
ITD = 255*10^-6;
% Parameters
f = 500;            % low frequency for MSO neurons.
fs = 20000;         % sampling rate
endTime = 0.01;      %end time of the sound wave in seconds
t = 0:1/fs:endTime;    % time axis

% Half-wave rectified inputs. Current assumption is that left input arrives
% first.
L = max(0, sin(2*pi*f*t));
R = max(0, sin(2*pi*f*(t - ITD)));

plot(t, L);
hold on;
plot(t, R);
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
title('Left vs Right Ear Sinusoids with 255 µs ITD');
legend('Left Ear','Right Ear');


end
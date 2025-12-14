function [] = jeffressModel()
% we use 8 neurons for simplicity purposes.
% negative ITD means input reaches right ear first
neuron_itds = [-255, -200, -100, -50, 0, 245, 250, 255];
left_delay = zeros(1, 8);
right_delay = zeros(1, 8);

% this is the ITD of the incoming sound. The required range for our model
% is [-255.102, +255.102 μs]

ITD = -255*10^-6;

% Parameters
f = 500;            % low frequency for MSO neurons.
fs = 20000;         % sampling rate
endTime = 0.01;      %end time of the sound wave in seconds
t = 0:1/fs:endTime;    % time axis

% Half-wave rectified inputs. Current assumption is that left input arrives
% first.
L = max(0, sin(2*pi*f*t));
R = max(0, sin(2*pi*f*(t - ITD)));

% we assume that either left_delay or right_delay for a given neuron is 0.
% and these delays are hard coded.
for i = 1:1:8
    left_delay(i) = max(0, neuron_itds(i));
    right_delay(i) = max(0, -neuron_itds(i));
end

y_vals = zeros(1, 8);
% now we simulate each neuron for the given input


function [y] = delaySignal(x, delay_sec, fs)
    % fs is 20,000 which means we divide 1s into 20,000 
    % intervals, each being equal to around 50μs. We multiply 20,000 by
    % delay_us seconds to figure out how many intervals does the delay
    % correspond to. That will be our delay starting index.
    delay_us = delay_sec * 1e-6;     % convert µs → seconds
    idx = round(delay_us * fs);

    y = zeros(size(x));
    if idx < length(x)
        y(idx+1:end) = x(1:end-idx);
    end
end


tau = 0.002;          % 2 ms membrane time constant
Vth = 0.5;            % spike threshold
dt = 1/fs;

V = 0;
spikeCount = 0;
spikeVals = zeros(1, 8);
for i = 1:1:8
    L_del = delaySignal(L, left_delay(i), fs);
    R_del = delaySignal(R, right_delay(i), fs);
    V = zeros(1, length(t));
    %hold on;
    for n = 1:length(t)-1
        gain = 1300; 
        I = gain*L_del(n) * R_del(n);
        V(n+1) = V(n) + dt * (-V(n)/tau + I);
        if V(n+1) >= Vth
            spikeCount = spikeCount + 1;
            V(n+1) = 0;   % reset
        end
    
    end
    spikeVals(i) = spikeCount;
    spikeCount= 0;
    plot(V);
end
%hold off;
%legend('-255μs', '-200μs', '-100μs', '-50μs', '0μs', '245μs', '250μs', '255μs')
%xlabel("Time (1 unit = 50μs)")
%ylabel("Neuron Activity (Membrane Potential)")
%title("LIF activity for neurons with different ITDs (no thresholding)")
%hold off;
x_vals = 1:1:8;
disp(spikeVals);
%plot(x_vals, spikeVals);
%xlabel("Neuron Index")
%ylabel("Number of spikes")
%title("Number of spikes for each MSO neuron (ITD = 0μs)")

% convert source azimuth angle to ITD
radius = 0.0875; % in metres
speed = 343; % in m/s

k = 2;   % number of neurons to display

[~, order] = sort(spikeVals, 'descend');
topIdx = order(1:k);
topITDs = neuron_itds(topIdx);

sumItds = sum(topITDs);

theta_list = zeros(1, k);
actual_pos = asin((speed*ITD)/radius);

for j = 1:k
    ITD_sec = topITDs(j) * 1e-6;
    arg = (speed * ITD_sec) / radius;
    arg = max(-1, min(1, arg));   % safety clamp
    theta_list(j) = asin(arg);
end


figure;
theta_plot = linspace(0, 2*pi, 200);
plot(sin(theta_plot), cos(theta_plot), 'k'); 
hold on;
axis equal;
grid on;

for j = 1:k
    sz = 5 + 2*spikeVals(topIdx(j));
    h_red =  plot(sin(theta_list(j)), cos(theta_list(j)), ...
         'ro', 'MarkerSize', sz, 'LineWidth', 2);
end
h_blue = plot(sin(actual_pos), cos(actual_pos), ...
              'bo', 'MarkerSize', 10, 'LineWidth', 2);

legend([h_red, h_blue], ...
       {'Top MSO neurons location prediction', 'Actual Sound Position'}, ...
       'Location', 'best');

title('Population Coding of Sound Azimuth');

end

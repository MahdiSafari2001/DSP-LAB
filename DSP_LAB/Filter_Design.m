%%
fs = 10000;  % Sampling frequency (Hz)

% Normalize frequency
fn = fs / 2;

% Butterworth filter design (1st order)
butter_cutoff = 100;  % Hz
[b_butter, a_butter] = butter(1, butter_cutoff / fn, 'low');

% Chebyshev filter design (1st order)
cheby_cutoff = 100;  % Hz
cheby_ripple = 1;  % dB
[b_cheby, a_cheby] = cheby1(1, cheby_ripple, cheby_cutoff / fn, 'low');
% Frequency response for Butterworth filter
[H_butter, w] = freqz(b_butter, a_butter, 'whole', fs);

% Frequency response for Chebyshev filter
[H_cheby, w] = freqz(b_cheby, a_cheby, 'whole', fs);

% Chebyshev + Butterworth 
butter_cheby = H_butter .* H_cheby;

% All-pass filter design (1st order)
% All-pass filter can be designed with a delay parameter.
ap_delay = 0.1;  % Seconds
theta = ap_delay * fs;
b_ap = [theta, -1];
a_ap = [1, -theta];

% Frequency response for All-pass filter
[H_ap, w] = freqz(b_ap, a_ap, 'whole', fs);

% Final frequency response
final_output = H_ap .* butter_cheby;
%%
% Plot frequency responses
figure;
subplot(2, 2, 1);
plot(w, 20*log10(abs(butter_cheby)));
title('Butterworth + Chebyshev (Amplitude)');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;

subplot(2, 2, 2);
plot(w,angle(H_cheby));
title('Butterworth + Chebyshev (Phase)');
xlabel('Frequency (Hz)');
ylabel('Phase (Radian)');
grid on;

subplot(2, 2, 3);
plot(w, 20*log10(abs(butter_cheby)));
title('Final Frequency Response (Amplitude)');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;

subplot(2, 2, 4);
plot(w,angle(final_output));
title('Final Frequency Response (Phase)');
xlabel('Frequency (Hz)');
ylabel('Phase (Radian)');
grid on;



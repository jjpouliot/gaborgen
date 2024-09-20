[amp, freqs, fftcomp] = freqtag_FFT3D(EEG.data, 500);



figure(5), plot(freqs(2:101), amp(:, 2:101))   %plot the all sensors 

ax = gca;         %editing the plot
ax.FontSize = 18; %set the font size
ax.Box = 'off';   %remove the box around the plot
xlabel('Frequency (Hz)'), ylabel('Amplitude (μV)') %label the axis


figure(6), plot(freqs(2:101), amp([31], 2:101))   %plot the all sensors 

ax = gca;         %editing the plot
ax.FontSize = 18; %set the font size
ax.Box = 'off';   %remove the box around the plot
xlabel('Frequency (Hz)'), ylabel('Amplitude (μV)') %label the axis

% after sliding window
[amp, freqs, fftcomp] = freqtag_FFT3D(winmat3d15Hz, resampledRateHz);

figure(7), bar(freqs(2:20), amp(:, 2:20))   %plot the all sensors 

ax = gca;         %editing the plot
ax.FontSize = 18; %set the font size
ax.Box = 'off';   %remove the box around the plot
xlabel('Frequency (Hz)'), ylabel('Amplitude (μV)') %label the axis


figure(8), bar(freqs(2:20), amp([31], 2:20))   %plot the all sensors 

ax = gca;         %editing the plot
ax.FontSize = 18; %set the font size
ax.Box = 'off';   %remove the box around the plot
xlabel('Frequency (Hz)'), ylabel('Amplitude (μV)') %label the axis
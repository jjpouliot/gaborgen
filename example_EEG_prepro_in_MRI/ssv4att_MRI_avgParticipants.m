partVec = [101:102 104:107];
avgChans = [9 10 20];

%% read and aggregate participant data

for partI = 1:length(partVec)
    if partI == 1
        dummyFFT857 = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_11_857_abs.txt']);
        dummyFFT15 = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_11_15_abs.txt']);
        dummyERP857 = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erp_11_857_cue.txt']);
        dummyERP15 = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erp_11_15_cue.txt']);
        dummyHilbert = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_hilbert_11_857.txt']);
        dummyPowST = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_fftST_11_857_abs.txt']);

        bigMatFFTabs857 = NaN(size(dummyFFT857,1),size(dummyFFT857,2),4,length(partVec));
        bigMatFFTabs15 = NaN(size(dummyFFT15,1),size(dummyFFT15,2),4,length(partVec));
        bigMatFFTblcorr857 = NaN(size(dummyFFT857,1),size(dummyFFT857,2),4,length(partVec));
        bigMatFFTblcorr15 = NaN(size(dummyFFT15,1),size(dummyFFT15,2),4,length(partVec));
        bigMatERP857 = NaN(size(dummyERP857,1),size(dummyERP857,2),4,length(partVec));
        bigMatERP15 = NaN(size(dummyERP15,1),size(dummyERP15,2),4,length(partVec));
        bigMatHilbert = NaN(size(dummyHilbert,1),size(dummyHilbert,2),8,length(partVec));
        bigMatPowST = struct();
        
        trialMat = NaN(length(partVec),4);
    end

    bigMatFFTabs857(:,:,1,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_11_857_abs.txt']); 
    bigMatFFTabs857(:,:,2,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_21_857_abs.txt']);   
    bigMatFFTabs857(:,:,3,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_22_857_abs.txt']);   
    bigMatFFTabs857(:,:,4,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_12_857_abs.txt']);

    bigMatFFTabs15(:,:,1,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_12_15_abs.txt']); 
    bigMatFFTabs15(:,:,2,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_22_15_abs.txt']);   
    bigMatFFTabs15(:,:,3,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_21_15_abs.txt']);   
    bigMatFFTabs15(:,:,4,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_11_15_abs.txt']);

    bigMatFFTblcorr857(:,:,1,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_11_857_blcorr.txt']); 
    bigMatFFTblcorr857(:,:,2,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_21_857_blcorr.txt']);   
    bigMatFFTblcorr857(:,:,3,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_22_857_blcorr.txt']);   
    bigMatFFTblcorr857(:,:,4,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_12_857_blcorr.txt']);

    bigMatFFTblcorr15(:,:,1,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_12_15_blcorr.txt']); 
    bigMatFFTblcorr15(:,:,2,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_22_15_blcorr.txt']);   
    bigMatFFTblcorr15(:,:,3,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_21_15_blcorr.txt']);   
    bigMatFFTblcorr15(:,:,4,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erpfft_11_15_blcorr.txt']);

    bigMatERP857(:,:,1,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erp_11_857_cue.txt']); 
    bigMatERP857(:,:,2,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erp_21_857_cue.txt']);   
    bigMatERP857(:,:,3,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erp_22_857_cue.txt']);   
    bigMatERP857(:,:,4,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erp_12_857_cue.txt']);

    bigMatERP15(:,:,1,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erp_12_15_cue.txt']); 
    bigMatERP15(:,:,2,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erp_22_15_cue.txt']);   
    bigMatERP15(:,:,3,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erp_21_15_cue.txt']);   
    bigMatERP15(:,:,4,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_erp_11_15_cue.txt']);

    bigMatHilbert(:,:,1,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_hilbert_11_857.txt']); 
    bigMatHilbert(:,:,2,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_hilbert_21_857.txt']);   
    bigMatHilbert(:,:,3,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_hilbert_22_857.txt']);   
    bigMatHilbert(:,:,4,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_hilbert_12_857.txt']);
    bigMatHilbert(:,:,5,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_hilbert_12_15.txt']); 
    bigMatHilbert(:,:,6,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_hilbert_22_15.txt']);   
    bigMatHilbert(:,:,7,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_hilbert_21_15.txt']);   
    bigMatHilbert(:,:,8,partI) = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_hilbert_11_15.txt']);

    bigMatPowST.attCol857{partI} = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_fftST_11_857_abs.txt']);
    bigMatPowST.attGray857{partI} = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_fftST_21_857_abs.txt']);
    bigMatPowST.ignCol857{partI} = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_fftST_22_857_abs.txt']);
    bigMatPowST.ignGray857{partI} = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_fftST_12_857_abs.txt']);
    bigMatPowST.attCol15{partI} = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_fftST_12_15_abs.txt']);
    bigMatPowST.attGray15{partI} = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_fftST_22_15_abs.txt']);
    bigMatPowST.ignCol15{partI} = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_fftST_21_15_abs.txt']);
    bigMatPowST.ignGray15{partI} = readmatrix(['/Volumes/TOSHIBA_4TB/ssV4att_MRI/data/' int2str(partVec(partI)) '/EEG/' int2str(partVec(partI)) '_fftST_11_15_abs.txt']);

    trialMat(partI,1) = size(bigMatPowST.attCol857{partI},2);
    trialMat(partI,2) = size(bigMatPowST.attCol15{partI},2);
    trialMat(partI,3) = size(bigMatPowST.attGray857{partI},2);
    trialMat(partI,4) = size(bigMatPowST.attGray15{partI},2);
end



%% plot spectra (abs)
figure;

plotBins = 12;
freq857 = 0:(60/7/4):((plotBins-1)*(60/7/4));
freq15 = 0:(15/4):((plotBins-1)*(15/4));

subplot(2,2,1);
hold on;
plot(freq857, squeeze(mean(mean(bigMatFFTabs857(avgChans,1:plotBins,1,:),4),1)), '-r');
plot(freq857, squeeze(mean(mean(bigMatFFTabs857(avgChans,1:plotBins,3,:),4),1)), '--r');
hold off;
title('Color at 8.57 Hz');
legend('att','ign');
ylim([0 0.04]);

subplot(2,2,2);
hold on;
plot(freq857, squeeze(mean(mean(bigMatFFTabs857(avgChans,1:plotBins,2,:),4),1)), '-b');
plot(freq857, squeeze(mean(mean(bigMatFFTabs857(avgChans,1:plotBins,4,:),4),1)), '--b');
hold off;
title('Gray at 8.57 Hz');
legend('att','ign');
ylim([0 0.04]);

subplot(2,2,3);
hold on;
plot(freq15, squeeze(mean(mean(bigMatFFTabs15(avgChans,1:plotBins,1,:),4),1)), '-r');
plot(freq15, squeeze(mean(mean(bigMatFFTabs15(avgChans,1:plotBins,3,:),4),1)), '--r');
hold off;
title('Color at 15 Hz');
legend('att','ign');
ylim([0 0.02]);

subplot(2,2,4);
hold on;
plot(freq15, squeeze(mean(mean(bigMatFFTabs15(avgChans,1:plotBins,2,:),4),1)), '-b');
plot(freq15, squeeze(mean(mean(bigMatFFTabs15(avgChans,1:plotBins,4,:),4),1)), '--b');
hold off;
title('Gray at 15 Hz');
legend('att','ign');
ylim([0 0.02]);


%% plot spectra topos (abs)
topovec_att_col_857 = mean(bigMatFFTabs857(:,5,1,:),4);
topovec_ign_col_857 = mean(bigMatFFTabs857(:,5,3,:),4);
topovec_diff_col_857 = topovec_att_col_857 - topovec_ign_col_857;

topovec_att_gray_857 = mean(bigMatFFTabs857(:,5,2,:),4);
topovec_ign_gray_857 = mean(bigMatFFTabs857(:,5,4,:),4);
topovec_diff_gray_857 = topovec_att_gray_857 - topovec_ign_gray_857;

topovec_att_col_15 = mean(bigMatFFTabs15(:,5,1,:),4);
topovec_ign_col_15 = mean(bigMatFFTabs15(:,5,3,:),4);
topovec_diff_col_15 = topovec_att_col_15 - topovec_ign_col_15;

topovec_att_gray_15 = mean(bigMatFFTabs15(:,5,2,:),4);
topovec_ign_gray_15 = mean(bigMatFFTabs15(:,5,4,:),4);
topovec_diff_gray_15 = topovec_att_gray_15 - topovec_ign_gray_15;

max857 = max([topovec_att_col_857;topovec_ign_col_857;topovec_att_gray_857;topovec_ign_gray_857]);
max15 = max([topovec_att_col_15;topovec_ign_col_15;topovec_att_gray_15;topovec_ign_gray_15]);

chanlocs = open('/Volumes/TOSHIBA_4TB/ssV4att_MRI/MATLAB code/chanlocs_ssv4att_MRI.mat'); chanlocs = chanlocs.chanlocs;
eeglab;

figure;
subplot(3,4,1); topoplot(topovec_att_col_857,chanlocs); title('att-col-857'); caxis([-max857 max857]); cbar;
subplot(3,4,5); topoplot(topovec_ign_col_857,chanlocs); title('ign-col-857'); caxis([-max857 max857]); cbar; 
subplot(3,4,9); topoplot(topovec_diff_col_857,chanlocs); title('diff-col-857'); caxis([-max857/2 max857/2]); cbar;

subplot(3,4,2); topoplot(topovec_att_gray_857,chanlocs); title('att-gray-857'); caxis([-max857 max857]); cbar;
subplot(3,4,6); topoplot(topovec_ign_gray_857,chanlocs); title('ign-cogray-857'); caxis([-max857 max857]); cbar;
subplot(3,4,10); topoplot(topovec_diff_gray_857,chanlocs); title('diff-gray-857'); caxis([-max857/2 max857/2]); cbar;

subplot(3,4,3); topoplot(topovec_att_col_15,chanlocs); title('att-col-15'); caxis([-max15 max15]); cbar;
subplot(3,4,7); topoplot(topovec_ign_col_15,chanlocs); title('ign-col-15'); caxis([-max15 max15]); cbar;
subplot(3,4,11); topoplot(topovec_diff_col_15,chanlocs); title('diff-col-15'); caxis([-max15/2 max15/2]); cbar;

subplot(3,4,4); topoplot(topovec_att_gray_15,chanlocs); title('att-gray-15'); caxis([-max15 max15]); cbar;
subplot(3,4,8); topoplot(topovec_ign_gray_15,chanlocs); title('ign-gray-15'); caxis([-max15 max15]); cbar;
subplot(3,4,12); topoplot(topovec_diff_gray_15,chanlocs); title('diff-gray-15'); caxis([-max15/2 max15/2]); cbar;


%% plot spectra (blcorr)
figure;

plotBins = 12;
freq857 = 0:(60/7/4):((plotBins-1)*(60/7/4));
freq15 = 0:(15/4):((plotBins-1)*(15/4));

subplot(2,2,1);
hold on;
plot(freq857, squeeze(mean(mean(bigMatFFTblcorr857(avgChans,1:plotBins,1,:),4),1)), '-r');
plot(freq857, squeeze(mean(mean(bigMatFFTblcorr857(avgChans,1:plotBins,3,:),4),1)), '--r');
hold off;
title('Color at 8.57 Hz');
legend('att','ign');
ylim([0 3]);

subplot(2,2,2);
hold on;
plot(freq857, squeeze(mean(mean(bigMatFFTblcorr857(avgChans,1:plotBins,2,:),4),1)), '-b');
plot(freq857, squeeze(mean(mean(bigMatFFTblcorr857(avgChans,1:plotBins,4,:),4),1)), '--b');
hold off;
title('Gray at 8.57 Hz');
legend('att','ign');
ylim([0 3]);

subplot(2,2,3);
hold on;
plot(freq15, squeeze(mean(mean(bigMatFFTblcorr15(avgChans,1:plotBins,1,:),4),1)), '-r');
plot(freq15, squeeze(mean(mean(bigMatFFTblcorr15(avgChans,1:plotBins,3,:),4),1)), '--r');
hold off;
title('Color at 15 Hz');
legend('att','ign');
ylim([0 3]);

subplot(2,2,4);
hold on;
plot(freq15, squeeze(mean(mean(bigMatFFTblcorr15(avgChans,1:plotBins,2,:),4),1)), '-b');
plot(freq15, squeeze(mean(mean(bigMatFFTblcorr15(avgChans,1:plotBins,4,:),4),1)), '--b');
hold off;
title('Gray at 15 Hz');
legend('att','ign');
ylim([0 3]);


%% plot ERPs
figure;

time857 = 0:1000/480:(4*1000/60*7-1);
time15 = 0:1000/480:(4*1000/15-1);

subplot(2,2,1);
hold on;
plot(time857, squeeze(mean(mean(bigMatERP857(avgChans,:,1,:),4),1)), '-r');
plot(time857, squeeze(mean(mean(bigMatERP857(avgChans,:,3,:),4),1)), '--r');
hold off;
title('Color at 8.57 Hz');
legend('att','ign');
ylim([-0.04 0.04]);

subplot(2,2,2);
hold on;
plot(time857, squeeze(mean(mean(bigMatERP857(avgChans,:,2,:),4),1)), '-b');
plot(time857, squeeze(mean(mean(bigMatERP857(avgChans,:,4,:),4),1)), '--b');
hold off;
title('Gray at 8.57 Hz');
legend('att','ign');
ylim([-0.04 0.04]);

subplot(2,2,3);
hold on;
plot(time15, squeeze(mean(mean(bigMatERP15(avgChans,:,1,:),4),1)), '-r');
plot(time15, squeeze(mean(mean(bigMatERP15(avgChans,:,3,:),4),1)), '--r');
hold off;
title('Color at 15 Hz');
legend('att','ign');
ylim([-0.02 0.02]);

subplot(2,2,4);
hold on;
plot(time15, squeeze(mean(mean(bigMatERP15(avgChans,:,2,:),4),1)), '-b');
plot(time15, squeeze(mean(mean(bigMatERP15(avgChans,:,4,:),4),1)), '--b');
hold off;
title('Gray at 15 Hz');
legend('att','ign');
ylim([-0.02 0.02]);



%% plot single-subject ERPs
figure;

time857 = 0:1000/480:(4*1000/60*7-1);
time15 = 0:1000/480:(4*1000/15-1);

subplot(2,4,1);
plot(time857, squeeze(mean(bigMatERP857(avgChans,:,1,:),1)));
title('attend Color at 8.57 Hz');
ylim([-0.08 0.08]);

subplot(2,4,2);
plot(time857, squeeze(mean(bigMatERP857(avgChans,:,2,:),1)));
title('attend Gray at 8.57 Hz');
ylim([-0.08 0.08]);

subplot(2,4,3);
plot(time857, squeeze(mean(bigMatERP857(avgChans,:,3,:),1)));
title('ignore Color at 8.57 Hz');
ylim([-0.08 0.08]);

subplot(2,4,4);
plot(time857, squeeze(mean(bigMatERP857(avgChans,:,4,:),1)));
title('ignore Gray at 8.57 Hz');
ylim([-0.08 0.08]);

subplot(2,4,5);
plot(time15, squeeze(mean(bigMatERP15(avgChans,:,1,:),1)));
title('attend Color at 15 Hz');
ylim([-0.03 0.03]);

subplot(2,4,6);
plot(time15, squeeze(mean(bigMatERP15(avgChans,:,2,:),1)));
title('attend Gray at 15 Hz');
ylim([-0.03 0.03]);

subplot(2,4,7);
plot(time15, squeeze(mean(bigMatERP15(avgChans,:,3,:),1)));
title('ignore Color at 15 Hz');
ylim([-0.03 0.03]);

subplot(2,4,8);
plot(time15, squeeze(mean(bigMatERP15(avgChans,:,4,:),1)));
title('ignore Gray at 15 Hz');
ylim([-0.03 0.03]);



%% scatter plots (abs)
figure;

subplot(2,2,1);
plot(1:2, squeeze(mean(bigMatFFTabs857(avgChans,5,[1 3],:),1)),'-o');
set(gca,'XTick',1:2,'XTickLabel',{'att','ign'})
title('Color at 8.57 Hz');
xlim([0.5 2.5]);
ylim([0 0.06]);

subplot(2,2,2);
plot(1:2, squeeze(mean(bigMatFFTabs857(avgChans,5,[2 4],:),1)),'-o');
set(gca,'XTick',1:2,'XTickLabel',{'att','ign'})
title('Gray at 8.57 Hz');
xlim([0.5 2.5]);
ylim([0 0.06]);

subplot(2,2,3);
plot(1:2, squeeze(mean(bigMatFFTabs15(avgChans,5,[1 3],:),1)),'-o');
set(gca,'XTick',1:2,'XTickLabel',{'att','ign'})
title('Color at 15 Hz');
xlim([0.5 2.5]);
ylim([0 0.03]);
legend({'101','102','104','105','106','107'});

subplot(2,2,4);
plot(1:2, squeeze(mean(bigMatFFTabs15(avgChans,5,[2 4],:),1)),'-o');
set(gca,'XTick',1:2,'XTickLabel',{'att','ign'})
title('Gray at 15 Hz');
xlim([0.5 2.5]);
ylim([0 0.03]);



%% scatter plots (baseline-corrected)
figure;

plotBins = 12;
freq857 = 0:(60/7/4):((plotBins-1)*(60/7/4));
freq15 = 0:(15/4):((plotBins-1)*(15/4));

subplot(2,2,1);
plot(1:2, squeeze(mean(bigMatFFTblcorr857(avgChans,5,[1 3],:),1)),'-o');
set(gca,'XTick',1:2,'XTickLabel',{'att','ign'})
title('Color at 8.57 Hz');
xlim([0.5 2.5]);
%ylim([0 0.06]);

subplot(2,2,2);
plot(1:2, squeeze(mean(bigMatFFTblcorr857(avgChans,5,[2 4],:),1)),'-o');
set(gca,'XTick',1:2,'XTickLabel',{'att','ign'})
title('Gray at 8.57 Hz');
xlim([0.5 2.5]);
%ylim([0 0.06]);

subplot(2,2,3);
plot(1:2, squeeze(mean(bigMatFFTblcorr15(avgChans,5,[1 3],:),1)),'-o');
set(gca,'XTick',1:2,'XTickLabel',{'att','ign'})
title('Color at 15 Hz');
xlim([0.5 2.5]);
%ylim([0 0.03]);

subplot(2,2,4);
plot(1:2, squeeze(mean(bigMatFFTblcorr15(avgChans,5,[2 4],:),1)),'-o');
set(gca,'XTick',1:2,'XTickLabel',{'att','ign'})
title('Gray at 15 Hz');
xlim([0.5 2.5]);
%ylim([0 0.03]);



%% plot Hilbert
figure;

time = -1598:1000/480:4000;


subplot(2,2,1);
hold on;
plot(time, squeeze(mean(mean(bigMatHilbert(avgChans,:,1,:),4),1)), '-r');
plot(time, squeeze(mean(mean(bigMatHilbert(avgChans,:,3,:),4),1)), '--r');
hold off;
title('Color at 8.57 Hz');
legend('att','ign');
%ylim([-0.04 0.04]);

subplot(2,2,2);
hold on;
plot(time, squeeze(mean(mean(bigMatHilbert(avgChans,:,2,:),4),1)), '-b');
plot(time, squeeze(mean(mean(bigMatHilbert(avgChans,:,4,:),4),1)), '--b');
hold off;
title('Gray at 8.57 Hz');
legend('att','ign');
%ylim([-0.04 0.04]);

subplot(2,2,3);
hold on;
plot(time, squeeze(mean(mean(bigMatHilbert(avgChans,:,5,:),4),1)), '-r');
plot(time, squeeze(mean(mean(bigMatHilbert(avgChans,:,7,:),4),1)), '--r');
hold off;
title('Color at 15 Hz');
legend('att','ign');
%ylim([-0.02 0.02]);

subplot(2,2,4);
hold on;
plot(time, squeeze(mean(mean(bigMatHilbert(avgChans,:,4,:),4),1)), '-b');
plot(time, squeeze(mean(mean(bigMatHilbert(avgChans,:,8,:),4),1)), '--b');
hold off;
title('Gray at 15 Hz');
legend('att','ign');
%ylim([-0.02 0.02]);
clear;
close all;

load mig25.mat;
migdata=X;

mig_image=abs(fftshift(fft2(migdata)));
mig_image=mig_image/max(max(mig_image));
 
% figure;colormap(jet);pcolor((mig_image));shading flat;
% colorbar;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BiMig = mig_image > 0.2;
% figure;colormap(gray);pcolor((BiMig));shading flat;
% colorbar;

A = size(BiMig);

for x = 1:A(1)
    for y = 1:A(2)
        if BiMig (x,y) == 0
            BiMig2 (x,y) = -1i;
        else
            BiMig2 (x,y) = 1i;
        end
    end
end

figure;
plot(BiMig2,'*');
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);

k=100;
noi = zeros(A(1),A(2));

%% from 0.1 to 5 

allSNR1 = zeros(1,length(0.1:0.01:0.8));
allBER1 = zeros(1,length(0.1:0.01:0.8));
count = 1;

for H=0.1:0.01:0.8
    
    for c=1:k
        for b=1:(A(1))
            noise = (0.1+ 0.01*count)/sqrt(2) * (randn(1,(A(2))) + 1i*randn(1,(A(2))));
            noi (b,:) = noi(b,:)+noise;
        end
    end
    
    SignalE = sum(sum(abs(BiMig2).^2));
    AvgNoiseE = ((sum(sum((abs(noi).^2))))*2)/k;

    noi = noi;
    Variance = reshape(noi,1,[]);
    var(Variance);
    noiseMig = BiMig2+noi;
    
    
    % figure;
    % plot(noiseMig,'*');
    % xlim([-1.5 1.5]);
    % ylim([-1.5 1.5]);
    
    
    SNR = 10*log10(SignalE/AvgNoiseE);
    allSNR1(count) = SNR;
    
    for x=1:A(1)
        for y=1:A(2)
            if imag(noiseMig(x,y)) >= 0
                receivedMig (x,y) = 1;
            else
                receivedMig (x,y) = 0;
            end
        end
    end
    
    if ((count == 1))
        figure;colormap(gray);pcolor((receivedMig));shading flat;
        colorbar;
        figure;
        plot(noiseMig,'*');
        xlim([-1.5 1.5]);
        ylim([-1.5 1.5]);

    end

    if ((count == 60))
        figure;colormap(gray);pcolor((receivedMig));shading flat;
        colorbar;
        plot(noiseMig,'*');
        xlim([-1.5 1.5]);
        ylim([-1.5 1.5]);
    end
    
    % figure;colormap(gray);pcolor((receivedMig));shading flat;
    % colorbar;
    
    BE = 0;
    
    for x=1:A(1)
        for y=1:A(2)
            if (receivedMig (x,y) ~= BiMig (x,y))
                BE = BE + 1;
            end
        end
    end
    
    SentB = A(1)*A(2);
    
    BER = BE/SentB;
    allBER1 (count) = BER;

    count = count+1;

end

%% to here

figure
plot(allBER1,allSNR1,'lineWidth',2);
ylabel('SNR (dB)');
xlabel('BER');
title('SNR vs. BER Scheme A');












%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


BiMig = mig_image < 0.3;
% figure;colormap(gray);pcolor((BiMig));shading flat;
% colorbar;

A = size(BiMig);

for x = 1:A(1)
    for y = 1:A(2)
        if BiMig (x,y) == 0
            BiMig2 (x,y) = -1i;
        else
            BiMig2 (x,y) = 1i;
        end
    end
end

% figure;
% plot(BiMig2,'*');
% xlim([-1.5 1.5]);
% ylim([-1.5 1.5]);

k=100;
noi = zeros(A(1),A(2));

%% from 0.1 to 5 

allSNR2 = zeros(1,length(0.1:0.01:0.8));
allBER2 = zeros(1,length(0.1:0.01:0.8));
count = 1;

for H=0.1:0.01:0.8
    
    
    for c=1:k
        for b=1:(A(1))
            noise = (0.1+ 0.01*count)/sqrt(2) * (randn(1,(A(2))) + 1i*randn(1,(A(2))));
            noi (b,:) = noi(b,:)+noise;
        end
    end
    
    SignalE = sum(sum(abs(BiMig2).^2));
    AvgNoiseE = ((sum(sum((abs(noi).^2))))*2)/k;

    noi = noi;
    Variance = reshape(noi,1,[]);
    var(Variance);
    noiseMig = BiMig2+noi;
    
    % figure;
    % plot(noiseMig,'*');
    % xlim([-1.5 1.5]);
    % ylim([-1.5 1.5]);
    
    SNR = 10*log10(SignalE/AvgNoiseE);
    allSNR2(count) = SNR;
    
    for x=1:A(1)
        for y=1:A(2)
            if imag(noiseMig(x,y)) >= 0
                receivedMig (x,y) = 1;
            else
                receivedMig (x,y) = 0;
            end
        end
    end
    
    if ((count == 1))
        figure;colormap(gray);pcolor((receivedMig));shading flat;
        colorbar;
    end

    if ((count == 60))
        figure;colormap(gray);pcolor((receivedMig));shading flat;
        colorbar;
    end
    
    % figure;colormap(gray);pcolor((receivedMig));shading flat;
    % colorbar;
    
    BE = 0;
    
    for x=1:A(1)
        for y=1:A(2)
            if (receivedMig (x,y) ~= BiMig (x,y))
                BE = BE + 1;
            end
        end
    end
    
    SentB = A(1)*A(2);
    
    BER = BE/SentB;
    allBER2 (count) = BER;

    count = count+1;

end

%% to here

figure
plot(allBER2,allSNR2,'lineWidth',2);
ylabel('SNR (dB)');
xlabel('BER');
title('SNR vs. BER Scheme B');





clear;
close all;

load mig25.mat;
migdata=X;

mig_image=abs(fftshift(fft2(migdata)));
mig_image=mig_image/max(max(mig_image));
 
figure;colormap(jet);pcolor((mig_image));shading flat;
colorbar;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BiMig = mig_image > 0.2;
figure;colormap(gray);pcolor((BiMig));shading flat;
title('Scheme A Send');
colorbar;

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
title('Coding');
%BiMig2V = reshape(BiMig2,1,A(1)*A(2));
%scatterplot(BiMig2V);
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);

noi = zeros(A(1),A(2));

%% from 0.1 to 5 

allSNR1 = zeros(1,length(-15:1:15));
allBER1 = zeros(1,length(-15:1:15));
count = 1;

for H=-15:1:15
    
    SignalE = sum(sum(BiMig2.^2));


    txwithnoise1=awgn(BiMig2,H,'measured');
    size1 = size(txwithnoise1);

    txwithnoise2=reshape(txwithnoise1,1,size1(1)*size1(2));

    rec=sum(txwithnoise2);

    AvgNoiseE = ((sum(sum((abs(txwithnoise1).^2))))*2);

    %Variance = reshape(noi,1,[]);
    %var(Variance);
    noiseMig = txwithnoise1;
    
    
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
        figure;
        colormap(gray);pcolor((receivedMig));shading flat;
        title('Scheme A High Noise');
        colorbar;
        figure;
        plot(noiseMig,'*');
        title('Scheme A High Noise');
        xlim([-1.5 1.5]);
        ylim([-1.5 1.5]);

    end

    if ((count == 25))
        figure;
        colormap(gray);pcolor((receivedMig));shading flat;
        title('Scheme A Low Noise');
        colorbar;
        figure;
        plot(noiseMig,'*');
        title('Scheme A Low Noise');
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

EbN0startdb = -15;
EbN0enddb = 10;
count = 1;

for EbN0dB = EbN0startdb:EbN0enddb
    EsN0 = (10^(EbN0dB/10));
    total=0;
    for PHI=0:0.001:(pi-pi/2)
        total=total+exp((-EsN0*((sin(pi/2))^2))/((sin(PHI))^2));
    end
    SER=total*0.001/pi;
    BER(count)=SER/2;

    Ebno = 10^(EbN0dB/10);
    BER(count) = qfunc(sqrt(2*Ebno));

    count=count+1;
end
figure

%plot(allSNR1,allBER1,'lineWidth',2);
semilogy((-15:1:15),allBER1,'lineWidth',2);
xlabel('SNR (dB)');
ylabel('BER');
title('SNR vs. BER Scheme A');
grid on;
hold on;
semilogy(EbN0startdb:EbN0enddb,BER, 'LineWidth', 1)










%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


BiMig = mig_image < 0.3;
figure;colormap(gray);pcolor((BiMig));shading flat;
title('Scheme B Send');
colorbar;

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
title('Coding');
%BiMig2V = reshape(BiMig2,1,A(1)*A(2));
%scatterplot(BiMig2V);
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);

noi = zeros(A(1),A(2));

%% from 0.1 to 5 

allSNR1 = zeros(1,length(-15:1:15));
allBER1 = zeros(1,length(-15:1:15));
count = 1;

for H=-15:1:15
    
    SignalE = sum(sum(BiMig2.^2));


    txwithnoise1=awgn(BiMig2,H,'measured');
    size1 = size(txwithnoise1);

    txwithnoise2=reshape(txwithnoise1,1,size1(1)*size1(2));

    rec=sum(txwithnoise2);

    AvgNoiseE = ((sum(sum((abs(txwithnoise1).^2))))*2);

    %Variance = reshape(noi,1,[]);
    %var(Variance);
    noiseMig = txwithnoise1;
    
    
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
        figure;
        colormap(gray);pcolor((receivedMig));shading flat;
        title('Scheme B High Noise');
        colorbar;
        figure;
        plot(noiseMig,'*');
        title('Scheme B High Noise');
        xlim([-1.5 1.5]);
        ylim([-1.5 1.5]);

    end

    if ((count == 25))
        figure;
        colormap(gray);pcolor((receivedMig));shading flat;
        title('Scheme B Low Noise');
        colorbar;
        figure;
        plot(noiseMig,'*');
        title('Scheme B Low Noise');
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

EbN0startdb = -15;
EbN0enddb = 10;
count = 1;

for EbN0dB = EbN0startdb:EbN0enddb
    Ebno = 10^(EbN0dB/10);
    BER(count) = qfunc(sqrt(2*Ebno));
    count=count+1;
end
figure

%plot(allSNR1,allBER1,'lineWidth',2);
semilogy((-15:1:15),allBER1,'lineWidth',2);
xlabel('SNR (dB)');
ylabel('BER');
title('SNR vs. BER Scheme B');
grid on;
hold on;
semilogy(EbN0startdb:EbN0enddb,BER, 'LineWidth', 1)

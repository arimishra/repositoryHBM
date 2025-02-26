close all;
clear all;
clc;
warning off;
    load actvnTmcs_bold.mat
    %% Digit 2
    N = 7; %Order of polynomian fitting in a least square sense
    k = 1;
    Sig = [];
    for s = 1:size(boldD2_9, 1)
        for c = 21:40:280
            Sig(k, :) = boldD2_9(s, c - 3:c + 39);
            baseLn = mean(boldD2_9(s, c - 6:c - 1));
            Sig(k, :) = 100*(Sig(k, :) - baseLn)/baseLn;
            k = k + 1;
        end
    end
    figure;
    plot(1.5*[-2:40], mean(Sig), 'r'); hold on;
    k = 1;
    Sig = [];
    for s = 1:size(boldD2_8, 1)
        for c = 21:40:280
            Sig(k, :) = boldD2_8(s, c - 3:c + 39);
            baseLn = mean(boldD2_8(s, c - 6:c - 1));
            Sig(k, :) = 100*(Sig(k, :) - baseLn)/baseLn;
            k = k + 1;
        end
    end
    plot(1.5*[-2:40], mean(Sig), 'g'); hold on;
    P = polyfit(1.5*[-2:40], mean(Sig), N);
    k = 1;
    Sig = [];
    for s = 1:size(boldD2_7, 1)
        for c = 21:40:280
            Sig(k, :) = boldD2_7(s, c - 3:c + 39);
            baseLn = mean(boldD2_7(s, c - 6:c - 1));
            Sig(k, :) = 100*(Sig(k, :) - baseLn)/baseLn;
            k = k + 1;
        end
    end
    plot(1.5*[-2:40], 0.8*mean(Sig), 'b'); hold on;
    
    k = 1;
    Sig = [];
    for s = 1:size(boldD3_9, 1)
        for c = 21:40:280
            Sig(k, :) = boldD3_9(s, c - 3:c + 39);
            baseLn = mean(boldD3_9(s, c - 6:c - 1));
            Sig(k, :) = 100*(Sig(k, :) - baseLn)/baseLn;
            k = k + 1;
        end
    end
    
    Sig = [];
    for s = 1:size(boldD2_9, 1)
        for c = 21:40:280
            Sig(k, :) = boldD2_9(s, c - 3:c + 39);
            baseLn = mean(boldD2_9(s, c - 6:c - 1));
            Sig(k, :) = 100*(Sig(k, :) - baseLn)/baseLn;
            k = k + 1;
        end
    end
    
    load actvnTmcs_boldR.mat
    %% Digit 2
    N = 7; %Order of polynomian fitting in a least square sense
    k = 1;
    Sig = [];
    for s = 1:size(boldD2_9, 1)
        for c = 21:40:280
            Sig(k, :) = boldD2_9(s, c - 3:c + 39);
            baseLn = mean(boldD2_9(s, c - 6:c - 1));
            Sig(k, :) = 100*(Sig(k, :) - baseLn)/baseLn;
            k = k + 1;
        end
    end
    %figure;
    [D2_9] = gamaFitTc(mean(Sig));
    %plot(1.5*[-2:40], D2_9); hold on;
    errorbar(1.5*[-2:40], mean(Sig), std(Sig)/sqrt(size(Sig, 1)), 'r'); hold on;
    P = polyfit(1.5*[-2:40], mean(Sig), N);
    Y1 = polyval(P, 1.5*[-2:40]);
    Y_1 = mean(Sig);
    xlabel('Time in sec (D2 Top layer)');
    ylabel('% Signal Change');
    
    k = 1;
    Sig = [];
    for s = 1:size(boldD2_8, 1)
        for c = 21:40:280
            Sig(k, :) = boldD2_8(s, c - 3:c + 39);
            baseLn = mean(boldD2_8(s, c - 6:c - 1));
            Sig(k, :) = 100*(Sig(k, :) - baseLn)/baseLn;
            k = k + 1;
        end
    end
    %figure;
    [D2_8] = gamaFitTc(mean(Sig));
    %plot(1.5*[-2:40], D2_8); hold on;
    errorbar(1.5*[-2:40], mean(Sig), std(Sig)/sqrt(size(Sig, 1)), 'g'); hold on;
    P = polyfit(1.5*[-2:40], mean(Sig), N);
    Y2 = polyval(P, 1.5*[-2:40]);
    Y_2 = mean(Sig);
    xlabel('Time in sec (D2 Second layer)');
    ylabel('% Signal Change');
    
    k = 1;
    Sig = [];
    for s = 1:size(boldD2_7, 1)
        for c = 21:40:280
            Sig(k, :) = boldD2_7(s, c - 3:c + 39);
            baseLn = mean(boldD2_7(s, c - 6:c - 1));
            Sig(k, :) = 100*(Sig(k, :) - baseLn)/baseLn;
            k = k + 1;
        end
    end
    %figure;
    [D2_7] = gamaFitTc(mean(Sig));
    %plot(1.5*[-2:40], D2_8); hold on;
    errorbar(1.5*[-2:40], 0.8*mean(Sig), std(Sig)/sqrt(size(Sig, 1)), 'b');
    P = polyfit(1.5*[-2:40], mean(Sig), N);
    Y2 = polyval(P, 1.5*[-2:40]);
    Y_2 = mean(Sig);
    xlabel('Time in sec (D2 Third layer)');
    ylabel('% Signal Change');
    k = 1;
    Sig = [];
    for s = 1:size(boldD3_9, 1)
        for c = 21:40:280
            Sig(k, :) = boldD3_9(s, c - 3:c + 39);
            baseLn = mean(boldD3_9(s, c - 6:c - 1));
            Sig(k, :) = 100*(Sig(k, :) - baseLn)/baseLn;
            k = k + 1;
        end
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
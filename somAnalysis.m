clear all; 
close all;
warning off;
    %% digitPSF_S(Stimulus area 3b) digitPSF_R(Resting area 3b) area1PSF_R(Resting area 1);
    %% stimulus FWHM
    btStrapN = 500;
    %{
    load dgtPSF_S.mat
    for s = 1:4
        if s == 1
            digit_Up = [digit2T; digit3T]';
            digit_Lr = [digit2B; digit3B]';
        elseif s == 2
            %% Include major minor axis only
            digit_Up = [digit2T(:, 1:2); digit3T(:, 1:2)]';
            digit_Lr = [digit2B(:, 1:2); digit3B(:, 1:2)]';
        elseif s == 3
            %% Include major minor axis and Area enclosed at FWHM
            digit_Up = [[digit2T(:, 1:2) digit2T(:, 4)]; [digit3T(:, 1:2) digit3T(:, 4)]]';
            digit_Lr = [[digit2B(:, 1:2) digit2B(:, 4)]; [digit3B(:, 1:2) digit3B(:, 4)]]';
            %% Include major minor axis and Axis ratio at FWHM
        else
            digit_Up = [[digit2T(:, 1:2) digit2T(:, 3)]; [digit3T(:, 1:2) digit3T(:, 3)]]';
            digit_Lr = [[digit2B(:, 1:2) digit2B(:, 3)]; [digit3B(:, 1:2) digit3B(:, 3)]]';
        end
        
        l1 = length(digit_Up);
        l2 = length(digit_Lr);
        for c = 1:btStrapN
            bstS1 = [];
            bstS2 = [];
            for i = 1:l1
                bstS1 = [bstS1 randperm(l1, 1)];
            end
            for i = 1:l2
                bstS2 = [bstS2 randperm(l2, 1)];
            end
            dataC = [digit_Up(:, bstS1) digit_Lr(:, bstS2)];
            %% SOM runs
            theclass = ones(l1 + l2, 1);
            theclass(l1 + 1:end) = 0;
            net = selforgmap([1 2]);
            net.trainParam.epochs = 200;
            net = train(net, dataC);
            idxSOM = vec2ind(net(dataC))';
            k(1) = sum(find(idxSOM(l1:end) == 1))/l1;
            k(2) = sum(find(idxSOM(l1 + 1:end) == 2))/l2;
            if k(1) > k(2)
                keeps1 = find(idxSOM == 1);
                keeps2 = find(idxSOM == 2);
                idxSOM(keeps1) = 0;
                idxSOM(keeps2) = 1;
            else
                keeps1 = find(idxSOM == 2);
                keeps2 = find(idxSOM == 1);
                idxSOM(keeps1) = 0;
                idxSOM(keeps2) = 1;
            end
            stimSOM_acr(c, s) = 100*length(find(idxSOM == theclass))/length(theclass);
        end
    end
    %}
    %%{
    %% resting FWHM
    load dgtPSF_R.mat
    for s = 1:9
        if s == 1
            digit_Up = [digit2T; digit3T]';
            digit_Lr = [digit2B; digit3B]';
        elseif s == 2
            %% Include major minor axis only
            digit_Up = [digit2T(:, 1:2); digit3T(:, 1:2)]';
            digit_Lr = [digit2B(:, 1:2); digit3B(:, 1:2)]';
        elseif s == 3
            %% Include major minor axis and Area enclosed at FWHM
            digit_Up = [[digit2T(:, 1:2) digit2T(:, 4)]; [digit3T(:, 1:2) digit3T(:, 4)]]';
            digit_Lr = [[digit2B(:, 1:2) digit2B(:, 4)]; [digit3B(:, 1:2) digit3B(:, 4)]]';
            %% Include major minor axis and Axis ratio at FWHM
        elseif s == 4
            digit_Up = [[digit2T(:, 1:2) digit2T(:, 3)]; [digit3T(:, 1:2) digit3T(:, 3)]]';
            digit_Lr = [[digit2B(:, 1:2) digit2B(:, 3)]; [digit3B(:, 1:2) digit3B(:, 3)]]';
        elseif s == 5
            digit_Up = [digit2T; digit3T]';
            digit_Lr = [digit2L; digit3L]';
        elseif s == 6
            %% Include major minor axis only
            digit_Up = [digit2T(:, 1:2); digit3T(:, 1:2)]';
            digit_Lr = [digit2L(:, 1:2); digit3L(:, 1:2)]';
        elseif s == 7
            %% Include major minor axis and Area enclosed at FWHM
            digit_Up = [[digit2T(:, 1:2) digit2T(:, 4)]; [digit3T(:, 1:2) digit3T(:, 4)]]';
            digit_Lr = [[digit2L(:, 1:2) digit2L(:, 4)]; [digit3L(:, 1:2) digit3L(:, 4)]]';
            %% Include major minor axis and Axis ratio at FWHM
        elseif s == 8
            digit_Up = [[digit2T(:, 1:2) digit2T(:, 3)]; [digit3T(:, 1:2) digit3T(:, 3)]]';
            digit_Lr = [[digit2L(:, 1:2) digit2L(:, 3)]; [digit3L(:, 1:2) digit3L(:, 3)]]';
        else
            digit_Up = [[digit2B(:, 1:2) digit2B(:, 4)]; [digit3B(:, 1:2) digit3B(:, 4)]]';
            digit_Lr = [[digit2L(:, 1:2) digit2L(:, 4)]; [digit3L(:, 1:2) digit3L(:, 4)]]';
        end
        
        l1 = length(digit_Up);
        l2 = length(digit_Lr);
        for c = 1:btStrapN
            bstS1 = [];
            bstS2 = [];
            for i = 1:l1
                bstS1 = [bstS1 randperm(l1, 1)];
            end
            for i = 1:l2
                bstS2 = [bstS2 randperm(l2, 1)];
            end
            dataC = [digit_Up(:, bstS1) digit_Lr(:, bstS2)];
            %% SOM runs
            theclass = ones(l1 + l2, 1);
            theclass(l1 + 1:end) = 0;
            net = selforgmap([1 2]);
            net.trainParam.epochs = 200;
            net = train(net, dataC);
            idxSOM = vec2ind(net(dataC))';
            k(1) = sum(find(idxSOM(l1:end) == 1))/l1;
            k(2) = sum(find(idxSOM(l1 + 1:end) == 2))/l2;
            if k(1) > k(2)
                keeps1 = find(idxSOM == 1);
                keeps2 = find(idxSOM == 2);
                idxSOM(keeps1) = 0;
                idxSOM(keeps2) = 1;
            else
                keeps1 = find(idxSOM == 2);
                keeps2 = find(idxSOM == 1);
                idxSOM(keeps1) = 0;
                idxSOM(keeps2) = 1;
            end
            restSOM_acr(c, s) = 100*length(find(idxSOM == theclass))/length(theclass);
        end
    end
    %}
    %% resting state connectivity based
    load restCnvt.mat
    %%{
    for s = 1:4
        if s == 1
            digit_Up = ar3bT;
            digit_Lr = ar3bB;
        elseif s == 2
            digit_Up = ar3bT;
            digit_Lr = ar3bL;
        elseif s == 3
            digit_Up = ar3bB;
            digit_Lr = ar3bL;
        else
            digit_Up = ar3bT;
            digit_Lr = [ar3bB ar3bL];
        end
        
        l1 = length(digit_Up);
        l2 = length(digit_Lr);
        for c = 1:btStrapN
            bstS1 = [];
            bstS2 = [];
            for i = 1:l1
                bstS1 = [bstS1 randperm(l1, 1)];
            end
            for i = 1:l2
                bstS2 = [bstS2 randperm(l2, 1)];
            end
            dataC = [digit_Up(:, bstS1) digit_Lr(:, bstS2)];
            %% SOM runs
            theclass = ones(l1 + l2, 1);
            theclass(l1 + 1:end) = 0;
            net = selforgmap([1 2]);
            net.trainParam.epochs = 200;
            net = train(net, dataC);
            idxSOM = vec2ind(net(dataC))';
            k(1) = sum(find(idxSOM(l1:end) == 1))/l1;
            k(2) = sum(find(idxSOM(l1 + 1:end) == 2))/l2;
            if k(1) > k(2)
                keeps1 = find(idxSOM == 1);
                keeps2 = find(idxSOM == 2);
                idxSOM(keeps1) = 0;
                idxSOM(keeps2) = 1;
            else
                keeps1 = find(idxSOM == 2);
                keeps2 = find(idxSOM == 1);
                idxSOM(keeps1) = 0;
                idxSOM(keeps2) = 1;
            end
            cnvtSOM_acr(c, s) = 100*length(find(idxSOM == theclass))/length(theclass);
        end
    end
    %}
    %%{
    for s = 1:4
        if s == 1
            digit_Up = ar1T;
            digit_Lr = ar1B;
        elseif s == 2
            digit_Up = ar1T;
            digit_Lr = ar1L;
        elseif s == 3
            digit_Up = ar1B;
            digit_Lr = ar1L;
        else
            digit_Up = ar1T;
            digit_Lr = [ar1B ar1L];
        end
        
        l1 = length(digit_Up);
        l2 = length(digit_Lr);
        for c = 1:btStrapN
            bstS1 = [];
            bstS2 = [];
            for i = 1:l1
                bstS1 = [bstS1 randperm(l1, 1)];
            end
            for i = 1:l2
                bstS2 = [bstS2 randperm(l2, 1)];
            end
            dataC = [digit_Up(:, bstS1) digit_Lr(:, bstS2)];
            %% SOM runs
            theclass = ones(l1 + l2, 1);
            theclass(l1 + 1:end) = 0;
            net = selforgmap([1 2]);
            net.trainParam.epochs = 200;
            net = train(net, dataC);
            idxSOM = vec2ind(net(dataC))';
            k(1) = sum(find(idxSOM(l1:end) == 1))/l1;
            k(2) = sum(find(idxSOM(l1 + 1:end) == 2))/l2;
            if k(1) > k(2)
                keeps1 = find(idxSOM == 1);
                keeps2 = find(idxSOM == 2);
                idxSOM(keeps1) = 0;
                idxSOM(keeps2) = 1;
            else
                keeps1 = find(idxSOM == 2);
                keeps2 = find(idxSOM == 1);
                idxSOM(keeps1) = 0;
                idxSOM(keeps2) = 1;
            end
            cnvtSOM_acr(c, s) = 100*length(find(idxSOM == theclass))/length(theclass);
        end
    end
    %}
    %save somAccuracy.mat stimSOM_acr restSOM_acr cnvtSOM_acr3b_1 cnvtSOM_acr1_3b;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
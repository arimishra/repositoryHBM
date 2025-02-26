%% spm_FCCnvty.m V7
% Steady-state functional connectivity analysis toolbox for SPM8/12.  This
% function does all the actual work.
%% Author:
%% Arabinda Mishra
% $Id: spm_stdcnvty.m 13 2021-01-14 10:40:03Z mishra $
%% Edited on 11/25/2024
% We can generate a visual report for each analysis, including slices of
% the connectivity map overlaid on the first functional image, and a
% histogram of the connectivity values.

function spm_FCCnvty(varargin)
    warning off;
    %%      varargin arguments are
    %%      1  imgs    (Image list)
    %%      2  roi     (ROI file)
    %%      3  TR      (Repetition time)
    %%      4  LPF     (Lowpass cutoff)
    %%      5  cbr     (Confounds branch)
    %%      6  maptyp  (Map type - mean ROI or single voxels)
    %%      7  ftsave  (Save filtered data flag)
    %%      8  svbeta  (Save beta images flag)
    %%      9  mask    (Mask any region, OPTIONAL)
    %%      10  princ   (principal components)
    %%      11 onset_t (Onset time for block design)
    %%      12 bld_cbv (Select the contrast BOLD/CBV)
    %%      13 sliding window width (Select the window width)
    %%      14 Uncorrected p value threshold
  
    imgs           = varargin{1};
    roi            = varargin{2};
    TR             = varargin{3};
    LPF            = varargin{4};
    cbr            = varargin{5};
    maptyp         = varargin{6};
    ftsave         = varargin{7};
    svbeta         = varargin{8};
    mask           = varargin{9};
    princ          = varargin{10};
    onset_t        = varargin{11};
    bld_cbv        = varargin{12};
    slWin_width    = varargin{13};
    p_ValueTh      = varargin{14};
    %% percentage threshold can be changed by the eqn. below
    if TR == 0
        disp('....');
        disp('......');
        disp('Execution terminated ....');
        disp('Enter TR ....');
        return;
    end
    per_thrld = 0.2;
    snr_thrld = 0.2;
    stat_4DFile = 0;
    [p f e] = fileparts(imgs{1});
    if e(1:3) == '.gz';
        gunzip([p filesep f e]);
        imgs{1} = [p filesep f];
        [p f e] = fileparts(imgs{1});
    end
    Vout = spm_vol(imgs{1});
    M = spm_read_vols(Vout);
    if sum(size(M) ~= 0) == 4
        [imgs] = fCnvty4_3Dconv(imgs{1});
        stat_4DFile = 1;
        zpFile = f;
    else
        file_name = [f e];
        path_name = [p filesep];
        n_digit = 0;
        while ~isempty(str2num(file_name(end - n_digit - 4)))
            n_digit = n_digit + 1;
        end
        c = str2num(file_name(end - n_digit - 3:end - 4)) - 1;
        while (exist([path_name file_name(1:end - n_digit - 4) num2str(c, ['%0' num2str(n_digit) '.0f']) '.nii']))
            c = c - 1;
        end
        init_img = c + 1;
        c = init_img;
        while (exist([path_name file_name(1:end - n_digit - 4) num2str(c, ['%0' num2str(n_digit) '.0f']) '.nii']))
            c = c + 1;
        end
        final_img = c - 1;
        brick_img = f(1:end - n_digit);
        for c = init_img:final_img
            imgs{c} = [path_name brick_img num2str(c, ['%0' num2str(n_digit) '.0f']) '.nii'];
        end
    end
    dispSlice = round(size(M, 3)/1.6);
    %% ROI process
    if ~isempty(roi)
        if ~isempty(roi{1})
            for i = 1:size(roi, 1)
                [p f e] = fileparts(roi{i});
                roi{i} = fullfile(p, [f e(1:4)]);
            end
            roi_file = f;
        end
    end
    if ~isempty(mask)
        if ~isempty(mask{1})
            for i = 1:size(mask, 1)
                [p f e] = fileparts(mask{i});
                mask{i} = fullfile(p, [f e(1:4)]);
            end
        end
    end
    if ~isempty(princ)
        if ~isempty(princ{1})
            for i = 1:size(princ, 1)
                [p f e] = fileparts(princ{i});
                princ{i} = fullfile(p, [f e(1:4)]);
            end
        end
    end
    %}
    if sum(onset_t) > 0
        disp('Analyzing STIMULUS data ....');
    else
        disp('Initiating rsfMRI connectivity analysis ....');
    end
    %% Image input
    if ~isempty(imgs{1}) == 1
        nimgs = length(imgs);
        [imgpath img1file img1ext] = fileparts(imgs{1});
        if ~isempty(roi)
            if ~isempty(roi{1})
                [imgpath_roi img1file_roi img1ext_roi] = fileparts(roi{1});
            end
        else
            if sum(onset_t) == 0
                disp('## WARNING ... No ROI/s defined ...');
            end
        end
        Vimgs = spm_vol(strvcat(imgs));
        %% Image masking - use external mask if supplied; if not, calculate a
        %% histogram-based threshold
        if ~isempty(mask)
            Vout = spm_vol(mask{1});
            Ymask = spm_read_vols(Vout);
            keeps = find(Ymask > 0);
            Ymask = (Ymask > 0);
        else
            Ymask = spm_read_vols(Vimgs(round(nimgs/2)));
            disp('Using histogram_based mask, Default...');
            imgthr = spm_antimode(Ymask);
            keeps = find(Ymask > imgthr/(100/per_thrld))';
            Ymask = (Ymask > imgthr/(100/per_thrld));
        end
        
        %% Confounds (global signal)
        if cbr.glob == 1
            glob = nan(nimgs, 1);
            for i = 1:nimgs
                Y = spm_read_vols(Vimgs(i));
                glob(i) = mean(Y(keeps));
            end
            disp('Removing global signal')
            zglob = zscore(glob);
        else
            disp('Not removing global signal')
            zglob = [];
        end
        %% High pass filter parameters discrete cosine components addition to regressors
        K.RT = TR;
        K.row = 1:nimgs;
        K.HParam = 100;
        K = spm_filter(K);
        hpf = K.X0;
        zhpf = zscore(hpf);
        zmot = [];
        %% Realignment parameters
        if ~isempty(cbr.mot)
            mot = load(cbr.mot{:});
            if length(mot) ~= nimgs
                disp('....');
                disp('......');
                disp('Execution terminated ....');
                disp('Check number of time points in motion parameter file ....');
                return;
            end
            zmot = zscore(mot);
        end
        %% Add principal component of the time course associated with muscle, csf etc. areas
        Yimgs = spm_read_vols(Vimgs);
        im_size = size(Yimgs);
        prin_c = [];
        if ~isempty(princ)
            for c = 1:length(princ)
                VM = spm_vol(princ{c});
                M = round(spm_read_vols(VM));
                kps = find(M > 0);
                Y_p = zeros(nimgs, length(kps));
                for i = 1:nimgs
                    DataPC = squeeze(Yimgs(:, :, :, i));
                    Y_p(i, :) = DataPC(kps);
                end
                Y_p = post_processMap(Y_p);
                [coefs, scr, lat, tsqrd, pecnt, mu] = pca(Y_p);
                comp_size = 1;
                while (sum(pecnt(1:comp_size)) < 70)
                    comp_size = comp_size + 1;
                end
                prin_c = scr(:, 1:comp_size)*coefs(:, 1:comp_size)' + repmat(mu, size(Y_p, 1), 1);
            end
            prin_c(isnan(prin_c)) = 0;
            %prin_c = zscore(prin_c);
            prin_c = zscore(mean(prin_c, 2));
        end
        %% Load data and reshape to 2D long form for easy processing
        disp('Loading data...')
        %% Get the seed ROI time series by finding non-zero voxels in the ROI mask
        if ~isempty(roi)
            if ~isempty(roi{1})
                if size(roi, 1) > 1
                    disp('Loading ROIs ....')
                    for i = 1:size(roi, 1)
                        Vroi = spm_vol(roi{i});
                        Yroi = round(spm_read_vols(Vroi));
                        kpROI{i} = find(Yroi > 0);
                        indroi = find(Yroi > 0);
                        if length(indroi) == 0
                            disp(['### WARNING ### Possibly ROI ' num2str(i) ' is not a proper mask']);
                        end
                        check_roi = find(indroi(1) == keeps);
                        if sum(check_roi) == 0
                            disp('### WARNING ### ROI selected outside the mask');
                        end
                    end
                else
                    Vroi = spm_vol(roi{:});
                    disp('Loading ROI data')
                    Yroi = round(spm_read_vols(Vroi));
                    indroi = find(Yroi > 0);
                    if length(indroi) > 0
                        Yseed = zeros(nimgs, length(indroi));
                    else
                        disp('### WARNING ### Possibly ROI is not a proper mask');
                    end
                    check_roi = find(indroi(1) == keeps);
                    if sum(check_roi) == 0
                        disp('### WARNING ### ROI selected outside the mask');
                    end
                end
            end
        end
        %% Include line below if you want to see maps with all frames (without excluding the beginning volumes) 
        Y_s = zeros(nimgs, length(keeps));
        for i = 1:nimgs
            Data_mat = squeeze(Yimgs(:, :, :, i));
            Y_s(i, :) = Data_mat(keeps);
        end
        Data_mat = zeros(im_size(1:3));
        %% Save raw snr map
        Y_s = post_processMap(Y_s);
        snr_raw = mean(Y_s)./std(Y_s);
        snr_raw = post_processMap(snr_raw);
        snrMask = (snr_raw > snr_thrld);
        Data_mat(keeps) = snrMask.*snr_raw;
        Vsnr = spm_vol(imgs{1});
        if sum(onset_t) > 0
            imgpath_roi = imgpath;
        end
        Vsnr.fname = [imgpath_roi filesep 'snrRaw.nii'];
        Vsnr = rmfield(Vsnr, 'pinfo');
        spm_write_vol(Vsnr, Data_mat);
        Xr = squeeze(Data_mat(:, :, dispSlice));
        
        if sum(onset_t) == 0
            roi_indx = check_roi(1);
        else
            roi_indx = round(length(keeps)/3);
        end
        
        %% Create confound removal design matrix Regresion starts
        disp('Confound removal...');
        dspTc = Y_s(:, roi_indx)';
        confX = [zmot prin_c zglob ones(nimgs, 1)];
        confX(isnan(confX)) = 0;
        b = lscov(confX, Y_s);
        Y_s = Y_s - confX(:, 1:end - 1)*b(1:end - 1, :);
        Y_s = post_processMap(Y_s);
        snr_reg = post_processMap(mean(Y_s)./std(Y_s));
        Data_mat(keeps) = snrMask.*snr_reg;
        Vsnr = spm_vol(imgs{1});
        if sum(onset_t) > 0
            imgpath_roi = imgpath;
        end
        Vsnr.fname = [imgpath_roi filesep 'snrReg.nii'];
        spm_write_vol(Vsnr, Data_mat);
      
        %% Block Design
        if sum(onset_t) > 0
            blockV = onset_t(3:4);
            length_intrv = onset_t(3:4);
            hrf = spm_hrf(TR);
            block = zeros(1, length_intrv(2));
            for s = 1:floor((nimgs + onset_t(1))/sum(blockV))
                block = [block ones(1, length_intrv(1)) zeros(1, length_intrv(2))];
            end
            if s*sum(blockV) + onset_t(2) ~= nimgs + onset_t(1)
                disp('....');
                disp('......');
                disp('Execution terminated ....');
                disp('CHECK PARADIGM ... ARE ON OFF BLOCKS SPECIFIED CORRECTLY? ...... ');
                return;
            end
            thisseed = block(onset_t(1) + 1:end)';
            if bld_cbv == 1
                k = 1;
                hrfCBV = [];
                for t = 1:TR:onset_t(3)*TR + 2*TR
                    hrfCBV(k) = -2.4486*(-0.184/1.5*exp(-t/1.5) + 0.330/4.5*exp(-t/4.5) + 0.670/13.5*exp(-t/13.5));
                    k = k + 1;
                end
                thisseed_hrf = conv(double(thisseed), hrfCBV);
                thisseed_hrf = thisseed_hrf(1:nimgs);
                hrf_wv = thisseed_hrf(onset_t(1) + 1:onset_t(1) + sum(onset_t(3:4)));
            else
                thisseed_hrf = conv(double(thisseed), hrf);
                thisseed_hrf = thisseed_hrf(1:nimgs);
                hrf_wv = thisseed_hrf(onset_t(1) + 1:onset_t(1) + sum(onset_t(3:4)));
            end
            pardmT = TR*(onset_t(2):sum(onset_t(3:4)):nimgs - onset_t(2)) - onset_t(1)*TR;
            seedGama = gamafittingHRF(hrf_wv')';
        end
        %% Additional specified confound time series
        %% RETROICOR Physiological signal correction
        if ~isempty(cbr.conf)
            [p f e] = fileparts(cbr.conf{:});
            scan_time = TR*nimgs; % Total Scan time
            win_wdth = 16;
            for c = 1:length(cbr.conf)
                if e(end) == 'g';
                    sampling_rate = 496;  %Sample rate for physiologic signal
                    no_data = sampling_rate*scan_time; %Total number of data point
                    fid = fopen(cbr.conf{c});
                    C_R = textscan(fid, '%n %n %n %n %n %n %n %n %n %n %n', 'commentStyle', '#');
                    phy_sig = [C_R{5} C_R{6} C_R{10}];
                    if no_data > size(phy_sig, 2)
                        disp('....');
                        disp('......');
                        disp('RETROICOR interupted: check physilogical signal file....');
                        disp('Recording might have stopped. Chech if recorded at 500Hz....');
                        return;
                    end
                    lastMark = max(find(squeeze(phy_sig(:, 3)) == 20));
                    phy_sig = phy_sig(lastMark - 2*no_data:2:end, :);
                    %phy_sig = phy_sig(length(phy_sig) - no_data + 1:end, :);
                    physiol_ecg = phy_sig(1:no_data, 1);
                    physiol_resp = phy_sig(1:no_data, 2);
                else
                    sampling_rate = 500;  %Sample rate for physiologic signal
                    no_data = sampling_rate*scan_time; %Total number of data point
                    fid = fopen(cbr.conf{c}, 'r');
                    phy_sig = fscanf(fid, '%f %f',[2 inf]);
                    fclose(fid);
                    if no_data > size(phy_sig, 2)
                        disp('....');
                        disp('......');
                        disp('RETROICOR interupted: check physilogical signal file....');
                        disp('Recording might have stopped. Chech if recorded at 500Hz....');
                        return;
                    end
                    physiol_ecg = phy_sig(1, 1:no_data)';
                    physiol_resp = phy_sig(2, 1:no_data)';
                end
                respCut_l = 0.1;
                respCut_h = 2;
                ecgCut_l = 2;
                ecgCut_h = 7;
                fs = scan_time/2 + 1/sampling_rate:1/sampling_rate:scan_time/2 + win_wdth;

                Rp = 3;   % Allowable passband ripple in dB
                Rs = 12;  % Required stopband attinuation in dB
                Wp = [2*respCut_l/sampling_rate 2*respCut_h/sampling_rate];
                Ws = [1.6*respCut_l/sampling_rate 2.4*respCut_h/sampling_rate];
                [n Wn] = cheb2ord(Wp, Ws, Rp, Rs);
                [b a] = cheby2(3, Rs, Wn);
                physiol_resp2 = filtfilt(b, a, physiol_resp);
                Wp = [2*ecgCut_l/sampling_rate 2*ecgCut_h/sampling_rate];
                Ws = [1.6*ecgCut_l/sampling_rate 2.4*ecgCut_h/sampling_rate];
                [n Wn] = cheb2ord(Wp, Ws, Rp, Rs);
                [b a] = cheby2(3, Rs, Wn);
                physiol_ecg2 = filtfilt(b, a, physiol_ecg);
                %% Spectral Analysis
                step_size = 10;
                frq_thrld = 10;
                Fs = 1/(2*step_size/1000);
                X = physiol_resp2(1:step_size:end);
                NFFT = 2^nextpow2(length(X));
                Y = fft(X, NFFT)/length(X);
                f = Fs/2*linspace(0, 1, NFFT/2 + 1);
                fs = round((2*frq_thrld/Fs)*length(f));
                P = abs(Y(1:NFFT/2 + 1));
                P = P/max(P(:));
                figure;
                subplot(2, 1, 1);
                plot(f(1:fs), P(1:fs));
                title(['Spectrum below Nyquist frequency ' num2str(Fs/2) 'Hz']);

                kps = find(P >= 0.3);
                kps1 = find(P == max(P));
                res_rblt = 100 - sum(abs(f(kps) - f(kps1)))*length(kps);
                xlabel('Frequency in Hz');
                ylabel('Signal Magnitude(Respiratory)');
                disp(['Peak respiratory signal power at : ' num2str(60*round(f(kps1)*100)/100) 'cyc/min']);
                if res_rblt < 50
                    disp('Respiratory Signal Unreliable to be Used ...');
                end

                X = physiol_ecg2(1:step_size:end);
                NFFT = 2^nextpow2(length(X));
                Y = fft(X, NFFT)/length(X);
                f = Fs/2*linspace(0, 1, NFFT/2 + 1);
                P = abs(Y(1:NFFT/2 + 1));
                P = P/max(P(:));
                subplot(2, 1, 2);
                plot(f(1:fs), P(1:fs));
                kps = find(P >= 0.3);
                kps1 = find(P == max(P));
                ecg_rblt = 100 - sum(abs(f(kps) - f(kps1)))*length(kps);
                xlabel('Frequency in Hz');
                ylabel('Signal Magnitude(Cardiac)');
                disp(['Peak cardiac signal power at : ' num2str(60*round(f(kps1)*100)/100) 'cyc/min']);
                ecg_Flag = 0;
                if ecg_rblt < 50
                    disp('Cardiac Signal Unreliable to be Used ... APPLYING RETROICOR ...');
                    disp('FOR Respiratory signal alone');
                else
                    disp('APPLYING RETROICOR ... for both cadiac and respiratory signal ...');
                    ecg_Flag = 1;
                end

                %% RetroIcor Physiological signal correction
                y = (physiol_resp2 > mean(physiol_resp2));
                c = 1;
                k = 1;
                cc = [];
                y(c) = ~y(c + 1);
                while (y(c)*y(c + 1) ~= 1 & c < length(y) - 1)
                    c = c + 1;
                    p = 0;
                    while (y(c)*y(c + 1) == 1 & c < length(y) - 1)
                        if p == 0
                            c1 = c;
                        end
                        cc(c) = c1;
                        c = c + 1;
                        p = p + 1;
                    end
                    if p ~= 0
                        cc(k) = round(c1 + p/2);
                        k = k + 1;
                    end
                end
                cc = cc(1:k - 1);
                y = y*0;
                y(cc) = 1;
                kps1 = find(y == 1);
                p = 1;
                phi_resp = [];
                k = 1;
                for c = 1.5*TR*sampling_rate + 1:sampling_rate*TR:no_data
                    while (c > kps1(p) & p < length(kps1));
                        p = p + 1;
                    end
                    if p > 1
                        phi_resp(k) = 2*pi*(c - kps1(p - 1))/(kps1(p) - kps1(p - 1));
                        k = k + 1;
                    end
                end

                %% spatial Phase correction (Respiratory)
                y_del_resp = [];
                for c = 1:size(Y_s, 2)
                    X = Y_s(:, c);
                    Xm = mean(X);
                    a1 = 0;
                    b1 = 0;
                    a2 = 0;
                    b2 = 0;
                    for n = 1:length(phi_resp)
                        a1 = a1 + 2*(X(n) - Xm)*cos(pi*phi_resp(n))/length(phi_resp);
                        b1 = b1 + 2*(X(n) - Xm)*sin(pi*phi_resp(n))/length(phi_resp);
                        a2 = a1 + 2*(X(n) - Xm)*cos(2*pi*phi_resp(n))/length(phi_resp);
                        b2 = b1 + 2*(X(n) - Xm)*sin(2*pi*phi_resp(n))/length(phi_resp);
                    end
                    for n = 1:length(phi_resp)
                        y_del(n + (size(Y_s, 1) - length(phi_resp))) = a1*cos(pi*phi_resp(n)) + a2*cos(2*pi*phi_resp(n)) + ...
                            b1*sin(pi*phi_resp(n)) + b2*sin(2*pi*phi_resp(n));
                    end
                    y_del_resp = [y_del_resp y_del'];
                end
                y_del_ecg = [];
                %% Cardiac signal correction
                if ecg_Flag == 1
                    y = (physiol_ecg2 > mean(physiol_ecg2));
                    c = 1;
                    k = 1;
                    cc = [];
                    y(c) = ~y(c + 1);
                    while (y(c)*y(c + 1) ~= 1 & c < length(y) - 1)
                        c = c + 1;
                        p = 0;
                        while (y(c)*y(c + 1) == 1 & c < length(y) - 1)
                            if p == 0
                                c1 = c;
                            end
                            cc(c) = c1;
                            c = c + 1;
                            p = p + 1;
                        end
                        if p ~= 0
                            cc(k) = round(c1 + p/2);
                            k = k + 1;
                        end
                    end
                    cc = cc(1:k - 1);
                    y = y*0;
                    y(cc) = 1;
                    kps1 = find(y == 1);
                    p = 1;
                    phi_ecg = [];
                    k = 1;
                    for c = 1.5*TR*sampling_rate + 1:sampling_rate*TR:no_data
                        while (c > kps1(p) & p < length(kps1));
                            p = p + 1;
                        end
                        if p > 1
                            phi_ecg(k) = 2*pi*(c - kps1(p - 1))/(kps1(p) - kps1(p - 1));
                            k = k + 1;
                        end
                    end
                    %% spatial Phase correction (Respiratory)
                    for c = 1:size(Y_s, 2)
                        X = Y_s(:, c);
                        Xm = mean(X);
                        a1 = 0;
                        b1 = 0;
                        a2 = 0;
                        b2 = 0;
                        for n = 1:length(phi_ecg)
                            a1 = a1 + 2*(X(n) - Xm)*cos(pi*phi_resp(n))/length(phi_ecg);
                            b1 = b1 + 2*(X(n) - Xm)*sin(pi*phi_resp(n))/length(phi_ecg);
                            a2 = a1 + 2*(X(n) - Xm)*cos(2*pi*phi_resp(n))/length(phi_ecg);
                            b2 = b1 + 2*(X(n) - Xm)*sin(2*pi*phi_resp(n))/length(phi_ecg);
                        end
                        for n = 1:length(phi_resp)
                            y_del(n + (size(Y_s, 1) - length(phi_resp))) = a1*cos(pi*phi_resp(n)) + a2*cos(2*pi*phi_resp(n)) + ...
                                b1*sin(pi*phi_resp(n)) + b2*sin(2*pi*phi_resp(n));
                        end
                        y_del_ecg = [y_del_ecg y_del'];
                    end
                    disp('Applying RETROICOR Cardiac signal correction switched off ...');
                    %% Y_s = Y_s - y_del_resp - y_del_ecg;
                    Y_s = Y_s - y_del_resp;
                else
                    Y_s = Y_s - y_del_resp;
                end
            end
        end
        dspTc1 = Y_s(:, roi_indx)';
        snr_ret = post_processMap(mean(Y_s)./std(Y_s));
        Data_mat(keeps) = snrMask.*snr_ret;
        Vsnr = spm_vol(imgs{1});
        if sum(onset_t) > 0
            imgpath_roi = imgpath;
        end
        Vsnr.fname = [imgpath_roi filesep 'snrRet.nii'];
        spm_write_vol(Vsnr, Data_mat);
        Xrr = squeeze(Data_mat(:, :, dispSlice));
        %% High pass/DC regression
        confX = [zhpf ones(nimgs, 1)];
        confX(isnan(confX)) = 0;
        b = lscov(confX, Y_s);
        Y_s = Y_s - confX(:, 1:end - 1)*b(1:end - 1, :);
        Y_s = post_processMap(Y_s);
        %% Band pass filteration of time series
        Rp = 3;   % Allowable passband ripple in dB
        Rs = 60;  % Required stopband attinuation in dB
        fs = 1/TR;
        if sum(onset_t) == 0
            if LPF > 0
                Wp = LPF/(fs/2);
                Ws = 1.2*Wp;
                if Ws > 1
                    disp('....');
                    disp('......');
                    disp('Execution terminated.... Upper cutoff more than nyquist frequency');
                    return;
                else
                    [n Wn] = cheb2ord(Wp, Ws, Rp, Rs);
                    [b a] = cheby2(n, Rs, Wn, 'low');
                    Y_s = filtfilt(b, a, Y_s);
                    Y_s = post_processMap(Y_s);
                end
            end
        else
            if fs > 0.5
                Wp = 0.25/(fs/2);
                Ws = 1.2*Wp;
                [n Wn] = cheb2ord(Wp, Ws, Rp, Rs);
                [b a] = cheby2(n, Rs, Wn);
                Y_s = filtfilt(b, a, Y_s);
                Y_s = post_processMap(Y_s);
            else
                disp('....');
                disp('Did not apply LPF, because Nyquist frequency is below 0.25Hz')
            end
        end
        Y_s = post_processMap(Y_s);
        %Y_s = mean(Y(:)) + detrend(Y_s);
        %% Pre-processed signal after LPF
        %%{
        Data_mat = zeros(im_size(1:3));
        for i = 1:nimgs
            Data_mat(keeps) = Y_s(i, :);
            Yimgs(:, :, :, i) = Data_mat;
        end
        
        dspTc2 = Y_s(:, roi_indx)';
        Data_mat = zeros(im_size(1:3));
        snr_lpf = mean(Y_s)./std(Y_s);
        Data_mat(keeps) = snrMask.*snr_lpf;
        Vsnr = spm_vol(imgs{1});
        if sum(onset_t) > 0
            imgpath_roi = imgpath;
        end
        Vsnr.fname = [imgpath_roi filesep 'snrPreprocessed.nii'];
        spm_write_vol(Vsnr, Data_mat);
        Xrrr = squeeze(Data_mat(:, :, dispSlice));
        if ~isempty(roi)
            if ~isempty(roi{1})
                if size(roi, 1) > 1
                    Yseed = [];
                    for c = 1:size(roi, 1)
                        indroi = kpROI{c};
                        xs = [];
                        for i = 1:nimgs
                            Data_mat = squeeze(Yimgs(:, :, :, i));
                            xs(i, :) = Data_mat(indroi);
                        end
                        ys{c} = xs;
                        if length(indroi) > 1
                            Yseed = [Yseed mean(xs')'];
                        else
                            Yseed = [Yseed xs];
                        end
                    end
                else
                    for i = 1:nimgs
                        Data_mat = squeeze(Yimgs(:, :, :, i));
                        Yseed(i, :) = Data_mat(indroi);
                    end
                end
            end
        end
        
        %% Write out filtered data with a new filename prefix 'prep_'
        %%{
        if ftsave == 1
            disp('Writing out pre-processed images');
            if stat_4DFile == 0
                [p f e] = fileparts(imgs{1});
                ofuncs = spm_select('FPList', [p filesep], ['^' brick_img '.*\.nii$']);
                matlabbatch{1}.spm.util.cat.name = [p filesep brick_img '4D.nii'];
                matlabbatch{1}.spm.util.cat.vols = cellstr(ofuncs);
                matlabbatch{1}.spm.util.cat.dtype = 4;
                spm_jobman('run', matlabbatch);
                %gzip(matlabbatch{1}.spm.util.cat.name);
                zF = matlabbatch{1}.spm.util.cat.name;
                delete([zF(1:end - 4) '.mat']);
                %delete(zF);
                for c = 1:nimgs
                    [p f e] = fileparts(imgs{c});
                    V = spm_vol(imgs{c});
                    V = rmfield(V, 'pinfo');
                    spm_write_vol(V, squeeze(Yimgs(:, :, :, c)));
                end
                ofuncs = spm_select('FPList', [p filesep], ['^' brick_img '.*\.nii$']);
                matlabbatch{1}.spm.util.cat.name = [p '/prep_' brick_img '4D.nii'];
                matlabbatch{1}.spm.util.cat.vols = cellstr(ofuncs);
                matlabbatch{1}.spm.util.cat.dtype = 4;
                spm_jobman('run', matlabbatch);
                %gzip(matlabbatch{1}.spm.util.cat.name);
                zF = matlabbatch{1}.spm.util.cat.name;
                delete([zF(1:end - 4) '.mat']);
                %delete(zF);
            else
                for c = 1:nimgs
                    [p f e] = fileparts(imgs{c});
                    V = spm_vol(imgs{c});
                    V = rmfield(V, 'pinfo');
                    spm_write_vol(V, squeeze(Yimgs(:, :, :, c)));
                end
                %% Writting zipped 4D nifti file
                disp('Writting Preprocessed 4D Image File ...');
                clear matlabbatch
                ofuncs = spm_select('FPList', [p filesep], '^rest.*\.nii$');
                matlabbatch{1}.spm.util.cat.name = [p '/prep_' zpFile '.nii'];
                matlabbatch{1}.spm.util.cat.vols = cellstr(ofuncs);
                matlabbatch{1}.spm.util.cat.dtype = 4;
                spm_jobman('run', matlabbatch);
                %gzip(matlabbatch{1}.spm.util.cat.name);
                zF = matlabbatch{1}.spm.util.cat.name;
                delete([zF(1:end - 4) '.mat']);
                %delete(zF);
            end
        else
            if stat_4DFile == 0
                for c = 1:nimgs
                    [p f e] = fileparts(imgs{c});
                    V = spm_vol(imgs{c});
                    V = rmfield(V, 'pinfo');
                    spm_write_vol(V, squeeze(Yimgs(:, :, :, c)));
                end
                [p f e] = fileparts(imgs{1});
                ofuncs = spm_select('FPList', [p filesep], ['^' brick_img '.*\.nii$']);
                matlabbatch{1}.spm.util.cat.name = [p filesep brick_img '4D.nii'];
                matlabbatch{1}.spm.util.cat.vols = cellstr(ofuncs);
                matlabbatch{1}.spm.util.cat.dtype = 4;
                spm_jobman('run', matlabbatch);
                %gzip(matlabbatch{1}.spm.util.cat.name);
                zF = matlabbatch{1}.spm.util.cat.name;
                delete([zF(1:end - 4) '.mat']);
                %delete(zF);
            else
                for c = 1:nimgs
                    [p f e] = fileparts(imgs{c});
                    V = spm_vol(imgs{c});
                    V = rmfield(V, 'pinfo');
                    spm_write_vol(V, squeeze(Yimgs(:, :, :, c)));
                end
                [p f e] = fileparts(imgs{1});
                ofuncs = spm_select('FPList', [p filesep], '^rest.*\.nii$');
            end
        end
        %}
        delete([imgpath_roi '/SPM.mat']);
        %% Can store the previous T map with as spmT_0001p.nii in case you want to compare two maps
        %{
        if exist([imgpath_roi '/spmT_0001.nii'])
            Vout = spm_vol([imgpath_roi '/spmT_0001.nii']);
            T = spm_read_vols(Vout);
            Vout.fname = [imgpath_roi filesep 'spmT_0001p.nii'];
            spm_write_vol(Vout, T);
        end
        %}
        
        clear zmot zglob zhpf conf prin_c M VM coefs data Y1;
        B = zeros(nimgs, max(size(keeps)));
        Z = zeros(nimgs, max(size(keeps)));
        pThreshold = p_ValueTh;
        Data_mat = zeros(im_size(1:3));
        if sum(onset_t) > 0
            %% spm12 GLM analysis
            %%{
            [p f e] = fileparts(imgs{1});
            clear matlabbatch
            matlabbatch{1}.spm.stats.fmri_spec.dir = {p};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(ofuncs);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond.name = 'Stimulus';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond.onset = pardmT;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond.duration = TR*blockV(1);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond.tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond.pmod = ...
                struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond.orth = 1;
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress = ...
                struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
            matlabbatch{1}.spm.stats.fmri_spec.fact = ...
                struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

            matlabbatch{2}.spm.stats.fmri_est.spmmat = {[p '/SPM.mat']};
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
            matlabbatch{3}.spm.stats.con.spmmat = matlabbatch{2}.spm.stats.fmri_est.spmmat;
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Pos Stimulus';
            if bld_cbv == 1
                matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [-1 0];
            else
                matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 0];
            end
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Neg Stimulus';
            if bld_cbv == 1
                matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 0];
            else
                matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 0];
            end
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.delete = 1;
            spm_jobman('run', matlabbatch);
            delete([p filesep zpFile '.mat']);
            %}
            crDat = zeros(sum(blockV), size(Y_s, 2));
            k1 = 1;
            for cc = onset_t(1) + 1:sum(blockV):nimgs - blockV(2)
                crDat = crDat + Y_s(cc:cc + sum(blockV) - 1, :);
                k1 = k1 + 1;
            end
            crDat = crDat/(k1 - 1);
            [cr p_val] = corr(crDat, hrf_wv);
            d1 = (p_val < pThreshold);
            d1 = d1.*snrMask';
            Data_mat(keeps) = d1.*cr;
            cr = Data_mat;
            Vout = spm_vol(imgs{1});
            Vout = rmfield(Vout, 'pinfo');
            if cbr.glob == 0
                Vout.fname = [imgpath filesep 'StimH.nii'];
                if exist([imgpath filesep 'StimH.nii']) == 2
                    disp('Over writting Stimulus driven activation maps in DATA FOLDER ...');
                end
            else
                Vout.fname = [imgpath filesep 'StimHG.nii'];
            end
            spm_write_vol(Vout, cr);
            
            X_hrf = [double(hrf_wv) ones(sum(blockV), 1)];
            b = lscov(X_hrf, crDat);
            Data_mat(keeps) = d1.*b(1, :)';
            B_r = Data_mat;
            Vout = spm_vol(imgs{1});
            Vout = rmfield(Vout, 'pinfo');
            if cbr.glob == 0
                Vout.fname = [imgpath filesep 'BstimH.nii'];
            else
                Vout.fname = [imgpath filesep 'BstimHG.nii'];
            end
            spm_write_vol(Vout, B_r);
            
            [cr p_val] = corr(crDat, seedGama);
            cr = d1.*cr;
            Zscr = atanh(cr)*sqrt(sum(blockV) - 1);
            Data_mat(keeps) = cr;
            cr = Data_mat;
            Vout = spm_vol(imgs{1});
            Vout = rmfield(Vout, 'pinfo');
            if cbr.glob == 0
                Vout.fname = [imgpath filesep 'Stim.nii'];
            else
                Vout.fname = [imgpath filesep 'StimG.nii'];
            end
            spm_write_vol(Vout, cr);
            Data_mat(keeps) = Zscr;
            Zscr = Data_mat;
            Vout = spm_vol(imgs{1});
            Vout = rmfield(Vout, 'pinfo');
            if cbr.glob == 0
                Vout.fname = [imgpath filesep 'Zstim.nii'];
            else
                Vout.fname = [imgpath filesep 'ZstimG.nii'];
            end
            spm_write_vol(Vout, Zscr);
            
            X_hrf = [double(seedGama) ones(sum(blockV), 1)];
            b = lscov(X_hrf, crDat);
            Data_mat(keeps) = d1.*b(1, :)';
            B_r = Data_mat;
            %B_r = B_r/max(B_r(:));
            Vout = spm_vol(imgs{1});
            Vout = rmfield(Vout, 'pinfo');
            if cbr.glob == 0
                Vout.fname = [imgpath filesep 'Bstim.nii'];
            else
                Vout.fname = [imgpath filesep 'BstimG.nii'];
            end
            spm_write_vol(Vout, B_r);
            
            [cr p_val] = corr(Y_s, thisseed_hrf);
            d1 = (p_val < pThreshold);
            p_val = d1.*p_val;
            Data_mat(keeps) = p_val;
            p_val = Data_mat;
            ks = find(p_val > 0);
            p_val(ks) = mafdr(p_val(ks), 'BHFDR', true);
            if cbr.glob == 0
                Vout.fname = [imgpath filesep 'pValue.nii'];
            else
                Vout.fname = [imgpath filesep 'pValueG.nii'];
            end
            spm_write_vol(Vout, p_val);
            Vout = spm_vol([imgpath filesep 'spmT_0001.nii']);
            P = (p_val < pThreshold).*spm_read_vols(Vout);
            spm_write_vol(Vout, P);
            Vout = spm_vol([imgpath filesep 'con_0001.nii']);
            P = (p_val < pThreshold).*spm_read_vols(Vout);
            spm_write_vol(Vout, P);
            Vout = spm_vol([imgpath filesep 'beta_0001.nii']);
            P = (p_val < pThreshold).*spm_read_vols(Vout);
            spm_write_vol(Vout, P);
            Vout = spm_vol([imgpath filesep 'spmT_0002.nii']);
            P = (p_val < pThreshold).*spm_read_vols(Vout);
            spm_write_vol(Vout, P);
            Vout = spm_vol([imgpath filesep 'con_0002.nii']);
            P = (p_val < pThreshold).*spm_read_vols(Vout);
            spm_write_vol(Vout, P);
            Vout = spm_vol([imgpath filesep 'beta_0002.nii']);
            P = (p_val < pThreshold).*spm_read_vols(Vout);
            spm_write_vol(Vout, P);
            
            Zscr = atanh(cr)*sqrt(nimgs - 3);
            Data_mat(keeps) = d1.*cr;
            cr = Data_mat;
            Vout = spm_vol(imgs{1});
            Vout = rmfield(Vout, 'pinfo');
            if cbr.glob == 0
                Vout.fname = [imgpath filesep 'StimT.nii'];
            else
                Vout.fname = [imgpath filesep 'StimTG.nii'];
            end
            spm_write_vol(Vout, cr);
            Data_mat(keeps) = d1.*Zscr;
            Zscr = Data_mat;
            Vout = spm_vol(imgs{1});
            Vout = rmfield(Vout, 'pinfo');
            if cbr.glob == 0
                Vout.fname = [imgpath filesep 'ZstimT.nii'];
            else
                Vout.fname = [imgpath filesep 'ZstimTG.nii'];
            end
            spm_write_vol(Vout, Zscr);
            
            X_hrf = [thisseed_hrf ones(nimgs, 1)];
            b = lscov(X_hrf, Y_s);
            Data_mat(keeps) = b(1, :)';
            B_r = Data_mat;
            Vout = spm_vol(imgs{1});
            Vout = rmfield(Vout, 'pinfo');
            if cbr.glob == 0
                Vout.fname = [imgpath filesep 'BstimT.nii'];
            else
                Vout.fname = [imgpath filesep 'BstimTG.nii'];
            end
            spm_write_vol(Vout, B_r);
        else
            if slWin_width == 0
                %% path to write maps in folder containing rois
                %imgpath = imgpath_roi
                %% Calculate mean seed time series if requested
                if size(roi, 1) > 1
                    crMat = ~eye(length(roi)).*corr(Yseed);
                    file_name = 'intrROIcor.mat';
                    corMat = zeros(length(ys) - 1);
                    for i = 1:length(ys)
                        for j = i + 1:length(ys)
                            cmt = corr([ys{i} ys{j}]);
                            corMat_r(i, j) = sum(sum(cmt(1:size(ys{i}, 2), size(ys{i}, 2) + 1:end)))/(size(ys{i}, 2)*size(ys{j}, 2));
                            corMat_val{i, j} = cmt(1:size(ys{i}, 2), size(ys{i}, 2) + 1:end);
                            corMat(i, j - 1) = crMat(i, j);
                            corZ(i, j) = atanh(corMat_r(i, j)) * sqrt(nimgs - 3);
                        end
                    end
                    corMat_r = corMat_r(1:i - 1, 2:i);
                    corZscr = corZ(1:i - 1, 2:i);
                    save([imgpath filesep file_name], 'Yseed', 'corMat', 'corMat_r', 'corMat_val', 'corZscr');
                    %save([imgpath filesep file_name], 'Yseed', 'corMat', 'corMat_r');
                else
                    if maptyp == 1 & size(Yseed, 2) > 1
                        Yseed = mean(Yseed')';
                    end
                    nseed = size(Yseed, 2);
                    Data_mat = zeros(im_size(1:3));
                    if  ~isempty(mask)
                        Vmask = spm_vol(mask{1});
                        Ymsk = spm_read_vols(Vmask);
                    else
                        Ymsk = ones(size(Data_mat));
                    end
                    for s = 1:nseed
                        thisseed = Yseed(:, s);
                        %[R p] = partialcorr(Y_ss, thisseed, confX);
                        [R p_val] = corr(Y_s, thisseed);
                        R = post_processMap(R);
                        p_val = post_processMap(p_val);
                        d1 = (p_val < pThreshold);
                        d1 = d1.*snrMask';
                        R = d1.*R;
                        Data_mat(keeps) = R;
                        r_value = Data_mat;
                        Zscr = atanh(R)*sqrt(nimgs - 3);
                        Data_mat(keeps) = Zscr;
                        Zscr = Ymask.*Data_mat;
                        Vout = spm_vol(imgs{1});
                        Vout = rmfield(Vout, 'pinfo');

                        if maptyp == 1 || nseed == 1
                            if cbr.glob == 0
                                Vout.fname = [imgpath filesep 'MnCnvt_' roi_file '.nii'];
                            else
                                Vout.fname = [imgpath filesep 'MnCnvtG_' roi_file '.nii'];
                            end
                            spm_write_vol(Vout, r_value);
                            if cbr.glob == 0
                                Vout.fname = [imgpath filesep 'ZMnCnvt_' roi_file '.nii'];
                            else
                                Vout.fname = [imgpath filesep 'ZMnCnvtG_' roi_file '.nii'];
                            end
                            spm_write_vol(Vout, Zscr);
                        else
                            if cbr.glob == 0
                                Vout.fname = [imgpath filesep 'ZCnvt(' roi_file ')' num2str(s, '%03.0f') '.nii'];
                            else
                                Vout.fname = [imgpath filesep 'ZCnvtG(' roi_file ')' num2str(s, '%03.0f') '.nii'];
                            end
                            spm_write_vol(Vout, Zscr);
                            if cbr.glob == 0
                                Vout.fname = [imgpath filesep 'cr_mat(' roi_file ')' num2str(s, '%03.0f') '.nii'];
                            else
                                Vout.fname = [imgpath filesep 'cr_amtG(' roi_file ')' num2str(s, '%03.0f') '.nii'];
                            end
                            spm_write_vol(Vout, r_value);
                        end
                        %% Bita values wrt seed voxel time course
                        if svbeta == 1
                            X = [double(thisseed) ones(nimgs, 1)];
                            b = lscov(X, Y_s);
                            B = post_processMap(b(1, :));
                            Data_mat(keeps) = d1.*B';
                            B_r = Ymsk.*Data_mat;
                            Vout = spm_vol(imgs{1});
                            Vout = rmfield(Vout, 'pinfo');
                            if maptyp == 1 || nseed == 1
                                if cbr.glob == 0
                                    Vout.fname = [imgpath filesep 'MnBCnvt_' roi_file '.nii'];
                                else
                                    Vout.fname = [imgpath filesep 'MnBCnvtG_' roi_file '.nii'];
                                end
                            else
                                if cbr.glob == 0
                                    Vout.fname = [imgpath filesep 'BCnvt(' roi_file ')' num2str(s, '%03.0f') '.nii'];
                                else
                                    Vout.fname = [imgpath filesep 'BCnvtG(' roi_file ')' num2str(s, '%03.0f') '.nii'];
                                end
                            end
                            spm_write_vol(Vout, B_r);
                        end
                    end

                    if maptyp == 0
                        brick_img = 'BCnvt';
                        ofuncs = spm_select('FPList', [imgpath filesep], ['^' brick_img '.*\.nii$']);
                        matlabbatch{1}.spm.util.cat.name = [imgpath '/Bcvt(' roi_file ')4D.nii'];
                        matlabbatch{1}.spm.util.cat.vols = cellstr(ofuncs);
                        matlabbatch{1}.spm.util.cat.dtype = 4;
                        spm_jobman('run', matlabbatch);
                        zF = matlabbatch{1}.spm.util.cat.name;
                        delete([zF(1:end - 4) '.mat']);

                        brick_img = 'ZCnvt';
                        ofuncs = spm_select('FPList', [imgpath filesep], ['^' brick_img '.*\.nii$']);
                        matlabbatch{1}.spm.util.cat.name = [imgpath '/Zcvt(' roi_file ')4D.nii'];
                        matlabbatch{1}.spm.util.cat.vols = cellstr(ofuncs);
                        matlabbatch{1}.spm.util.cat.dtype = 4;
                        spm_jobman('run', matlabbatch);
                        zF = matlabbatch{1}.spm.util.cat.name;
                        delete([zF(1:end - 4) '.mat']);

                        brick_img = 'cr_mat';
                        ofuncs = spm_select('FPList', [imgpath filesep], ['^' brick_img '.*\.nii$']);
                        matlabbatch{1}.spm.util.cat.name = [imgpath '/Crcvt(' roi_file ')4D.nii'];
                        matlabbatch{1}.spm.util.cat.vols = cellstr(ofuncs);
                        matlabbatch{1}.spm.util.cat.dtype = 4;
                        spm_jobman('run', matlabbatch);
                        zF = matlabbatch{1}.spm.util.cat.name;
                        delete([zF(1:end - 4) '.mat']);
                    end

                    for c = 1:nseed
                        delete([imgpath '/BCnvt(' roi_file ')' num2str(c, '%03.0f') '.nii']);
                        delete([imgpath '/ZCnvt(' roi_file ')' num2str(c, '%03.0f') '.nii']);
                        delete([imgpath '/cr_mat(' roi_file ')' num2str(c, '%03.0f') '.nii']);
                    end

                    if size(Yseed, 2) > 1
                        Yseed = mean(Yseed')';
                    end
                    [cr p_val] = corr(Y_s, Yseed);
                    p_val = post_processMap(p_val);
                    d1 = (p_val < pThreshold);
                    d1 = d1.*snrMask';
                    Data_mat(keeps) = d1.*p_val;
                    p_val = Data_mat;
                    ks = find(p_val > 0);
                    p_val(ks) = mafdr(p_val(ks),'BHFDR', true);
                    Vout = spm_vol(imgs{1});
                    Vout = rmfield(Vout, 'pinfo');
                    if cbr.glob == 0
                        Vout.fname = [imgpath filesep 'pValue.nii'];
                    else
                        Vout.fname = [imgpath filesep 'pValueG.nii'];
                    end
                    spm_write_vol(Vout, p_val);
                end
            else
                imgpath = [imgpath(1:end - 8) 'anat_img'];
                if size(Yseed, 2) > 1
                    Yseed = mean(Yseed')';
                end
                stepSize = 2;
                if exist([imgpath filesep 'tempDyn4D.nii']) == 2
                    disp('Over writting temporal connectivity Z maps ...');
                end
                k = 1;
                for c = 1:stepSize:nimgs - slWin_width
                    thisseed = Yseed(c:c + slWin_width - 1);
                    [R p_val] = corr(Y_s(c:c + slWin_width - 1, :), thisseed);
                    p_val = post_processMap(p_val);
                    d1 = (p_val < pThreshold);
                    d1 = d1.*snrMask';
                    R = d1.*R;
                    R = post_processMap(R);
                    Zscr = atanh(R)*sqrt(nimgs - 3);
                    Zscr = d1.*Zscr;
                    Data_mat(keeps) = R;
                    r_value = Data_mat;
                    Data_mat(keeps) = Zscr;
                    Zscr = Data_mat;
                    Vout = spm_vol(imgs{1});
                    Vout = rmfield(Vout, 'pinfo');
                    Vout.fname = [imgpath filesep 'TempDynmcs(' roi_file ')' num2str(k, '%03.0f') '.nii'];
                    spm_write_vol(Vout, r_value);
                    Vout.fname = [imgpath filesep 'ZTempDynmcs(' roi_file ')' num2str(k, '%03.0f') '.nii'];
                    spm_write_vol(Vout, Zscr);
                    k = k + 1;
                end

                brick_img = 'TempDynmcs';
                ofuncs = spm_select('FPList', [imgpath filesep], ['^' brick_img '.*\.nii$']);
                matlabbatch{1}.spm.util.cat.name = [imgpath '/tempDyn4D.nii'];
                matlabbatch{1}.spm.util.cat.vols = cellstr(ofuncs);
                matlabbatch{1}.spm.util.cat.dtype = 4;
                spm_jobman('run', matlabbatch);
                zF = matlabbatch{1}.spm.util.cat.name;
                delete([zF(1:end - 4) '.mat']);

                brick_img = 'ZTempDynmcs';
                ofuncs = spm_select('FPList', [imgpath filesep], ['^' brick_img '.*\.nii$']);
                matlabbatch{1}.spm.util.cat.name = [imgpath '/zTempDyn4D.nii'];
                matlabbatch{1}.spm.util.cat.vols = cellstr(ofuncs);
                matlabbatch{1}.spm.util.cat.dtype = 4;
                spm_jobman('run', matlabbatch);
                zF = matlabbatch{1}.spm.util.cat.name;
                delete([zF(1:end - 4) '.mat']);
            end
        end
    end
    
    %% Delete 3D volumes
    %%{
    if stat_4DFile == 1
        [p f e] = fileparts(imgs{1});
        for c = 1:nimgs
            delete([p '/rest' num2str(c, '%04.0f') '.nii']);
        end
    elseif exist([zpFile '.gz'], 'file')
        for c = init_img:final_img
            delete([path_name brick_img num2str(c, ['%0' num2str(n_digit) '.0f']) '.nii']);
        end
    end
    for c = 1:nimgs
        delete([imgpath '/Res_' num2str(c, '%04.0f') '.nii']);
    end
    %}

    %% Display the time course of four representative voxels.
    figure;
    subplot(2, 2, 1), plot(dspTc, 'r', 'LineWidth', 1), hold on;
    xlabel('Time Course: RAW DATA', 'fontname', 'arial', 'fontsize', 14);
    title('REPRESENTATIVE VOXEL TIME COURSE');
    set(gca,'FontSize', 14);

    subplot(2, 2, 2), plot(dspTc1, 'b', 'LineWidth', 1), hold on;
    xlabel('Time Course: PRE-PROCESSED', 'fontname', 'arial', 'fontsize', 14);
    set(gca,'FontSize', 14);

    subplot(2, 2, 3), plot(dspTc2, 'g', 'LineWidth', 1), hold on;
    xlabel('Time Course: LOW PASS FILTERED', 'fontname', 'arial', 'fontsize', 14);
    set(gca,'FontSize', 14);

    subplot(2, 2, 4), plot(dspTc, 'r', 'LineWidth', 1), hold on;
    subplot(2, 2, 4), plot(dspTc1, 'b', 'LineWidth', 1), hold on;
    subplot(2, 2, 4), plot(dspTc2, 'g', 'LineWidth', 1), hold on;
    xlabel('TIME COURSE OVERLAID', 'fontname', 'arial', 'fontsize', 14);
    set(gca,'FontSize', 14);

    figure;
    imagesc([Xr Xrr Xrrr]);
    title('SNR COMPARISON                             1.Raw                          2.Regressed                      3.LFP');
    set(gca,'FontSize', 14);
    axis off;
    colormap(jet);
    colorbar;
    saveas(gcf, [imgpath '/SNRcomp.png']);
    
    %% Signal magnitude vs. frequency analysis
    %{
    X = detrend(dspTc);
    Fs = 1/TR;
    frq_thrld = Fs/2;
    NFFT = 2^nextpow2(length(X));
    Y = fft(X, NFFT)/length(X);
    f = Fs/2*linspace(0, 1, NFFT/2 + 1);
    fs = round((2*frq_thrld/Fs)*length(f));
    P = abs(Y(1:NFFT/2 + 1));
    P = P/max(P(:));
    figure;
    plot(f(1:fs), P(1:fs), 'r', 'LineWidth', 1); hold on;
    
    X = detrend(dspTc1);
    NFFT = 2^nextpow2(length(X));
    Y = fft(X, NFFT)/length(X);
    f = Fs/2*linspace(0, 1, NFFT/2 + 1);
    fs = round((2*frq_thrld/Fs)*length(f));
    P = abs(Y(1:NFFT/2 + 1));
    P = P/max(P(:));
    plot(f(1:fs), P(1:fs), 'b', 'LineWidth', 1); hold on;
    
    X = detrend(dspTc2);
    NFFT = 2^nextpow2(length(X));
    Y = fft(X, NFFT)/length(X);
    f = Fs/2*linspace(0, 1, NFFT/2 + 1);
    fs = round((2*frq_thrld/Fs)*length(f));
    P = abs(Y(1:NFFT/2 + 1));
    P = P/max(P(:));
    plot(f(1:fs), P(1:fs), 'g', 'LineWidth', 1);
    title(['Spectrum below Nyquist frequency (red-Raw, black-regressed, green-LPF ...' num2str(Fs/2) 'Hz']);
    set(gca,'FontSize', 14);
    %}
    
    %% Signal power vs. frequency (Spectral Power)
    %%{
    P = pwelch(detrend(dspTc));
    P = P/max(P);
    sZ = welchInterp(length(P), 0.5/TR);
    figure;
    plot(0:sZ:0.5/TR, P, 'r', 'LineWidth', 2); hold on;
    P = pwelch(detrend(dspTc1));
    P = P/max(P);
    plot(0:sZ:0.5/TR, P, 'b', 'LineWidth', 2); hold on;
    P = pwelch(detrend(dspTc2));
    P = P/max(P);
    plot(0:sZ:0.5/TR, P, 'g', 'LineWidth', 2);
    title(['Normalized Power Spectrum red-Raw, black-regressed, green-LPF']);
    set(gca,'FontSize', 14);
    %}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
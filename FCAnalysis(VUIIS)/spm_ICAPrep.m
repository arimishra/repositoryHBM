%% spm_stdcvtyICA.m
% Pre-process the data for functional connectivity analysis etc. toolbox for SPM12.
% Authors:
%%     Arabinda Mishra
%% April 20, 2015 Edited 11/25/2024
% Updated Jan 14th, 2021
% We can generate a visual report for each analysis, including slices of
% the connectivity map overlaid on the first functional image, and a
% histogram of the connectivity values.

function spm_ICAPrep(varargin)
    warning off;
    %%  varargin arguments are
    %%   1  imgs    (Image list)
    %%   2  TR      (Repetition time)
    %%   3  LPF     (Lowpass cutoff)
    %%   4  cbr     (Confounds branch)
    %%   5  mask    (Mask any region preferaby for muscle mask, white or gray mask regression(OPTIONAL)
    %%   6 princ   (principal components)
    %%   7 onset time for stimulus data pre-processing
    
    imgs           = varargin{1};
    TR             = varargin{2};
    LPF            = varargin{3};
    cbr            = varargin{4};
    mask           = varargin{5};
    princ          = varargin{6};
    onset_t        = varargin{7};
    
    %% percentage threshold can be changed by the eqn. below
    per_thrld = 0.2;
    snr_thrld = 0.2;
    stat_4DFile = 0;
    clc;
    if TR == 0
        disp('....');
        disp('......');
        disp('Execution terminated ....');
        disp('Enter TR ....');
        return;
    end
    if ~isempty(mask)
        [p f e] = fileparts(mask{1});
        mask = fullfile(p, [f e(1:4)]);
    end
    [p f e] = fileparts(imgs{:});
    if e(1:3) == '.gz';
        gunzip([p '/' f e]);
        imgs{1} = [p '/' f];
        [p f e] = fileparts(imgs{1});
    end
    Vout = spm_vol(imgs{1});
    M = spm_read_vols(Vout);
    if sum(size(M) ~= 0) == 4
        %gzip(imgs{1});
        [imgs] = fCnvty4_3Dconv(imgs{1});
        stat_4DFile = 1;
        zpFile = f;
        %delete([p '/' f e]);
    else
        file_name = [f e];
        path_name = [p '/'];
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
    Vimgs = spm_vol(strvcat(imgs));
    nimgs = length(imgs);
    [imgpath img1file img1ext] = fileparts(imgs{1});
    disp('Analysis initiated ....');
    dispSlice = round(size(M, 3)/1.6);
    %% Image input
    if ~isempty(mask)
        Vout = spm_vol(mask);
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
    %% Confounds
    %% Global signal
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
    else
        prin_c = [];
    end
    %% Load data and reshape to 2D long form for easy processing
    disp('Loading data...')
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
    imgpath_roi = imgpath;
    Vsnr.fname = [imgpath_roi filesep 'snrRaw.nii'];
    Vsnr = rmfield(Vsnr, 'pinfo');
    spm_write_vol(Vsnr, Data_mat);
    Xr = squeeze(Data_mat(:, :, dispSlice));
    %% Input index (roi_indx) to any to seee the spectra at any location
    roi_indx = round(length(keeps)/3);

    %% Create confound removal design matrix Regresion starts
    disp('Confound removal...');
    dspTc = Y_s(:, roi_indx)';
    %confX = [zmot prin_c zhpf zglob ones(nimgs, 1)];
    confX = [zmot prin_c zglob ones(nimgs, 1)];
    confX(isnan(confX)) = 0;
    b = lscov(confX, Y_s);
    Y_s = Y_s - confX(:, 1:end - 1)*b(1:end - 1, :);
    Y_s = post_processMap(Y_s);
    snr_reg = post_processMap(mean(Y_s)./std(Y_s));
    Data_mat(keeps) = snrMask.*snr_reg;
    Vsnr = spm_vol(imgs{1});
    Vsnr.fname = [imgpath_roi filesep 'snrReg.nii'];
    Vsnr = rmfield(Vsnr, 'pinfo');
    spm_write_vol(Vsnr, Data_mat);
    
    if ~isempty(cbr.conf)
        [p f e] = fileparts(cbr.conf{:});
        scan_time = TR*nimgs;
        win_wdth = 16;
        for c = 1:length(cbr.conf)
            if e(end) == 'g';
                sampling_rate = 496;  %Sample rate for physiologic signal
                no_data = sampling_rate*scan_time; %Total number of data point
                fid = fopen(cbr.conf{c});
                C_R = textscan(fid, '%n %n %n %n %n %n %n %n %n %n %n','commentStyle', '#');
                phy_sig = [C_R{5} C_R{6}];
                if no_data > size(phy_sig, 2)
                    disp('....');
                    disp('......');
                    disp('RETROICOR interupted: check physilogical signal file....');
                    disp('Recording might have stopped. Chech if recorded at 500Hz....');
                    return;
                end
                phy_sig = phy_sig(length(phy_sig) - no_data + 1:end, :);
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
            title('Signal Magnitude VS sampling frequency (Resp/Cardiac)');

            kps = find(P >= 0.3);
            kps1 = find(P == max(P));
            res_rblt = 100 - sum(abs(f(kps) - f(kps1)))*length(kps);
            xlabel('Frequency in Hz');
            ylabel('Signal Magnitude(Respiratory)');
            disp(['Peak respiratory signal power at : ' num2str(60*round(f(kps1)*100)/100) 'cyc/min']);

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
            if ecg_rblt < 50 | round(f(kps1)*100)/100 > 3*(1/TR)
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
    %}
    dspTc1 = Y_s(:, roi_indx)';
    snr_ret = post_processMap(mean(Y_s)./std(Y_s));
    Data_mat(keeps) = snrMask.*snr_ret;
    Vsnr = spm_vol(imgs{1});
    Vsnr.fname = [imgpath_roi filesep 'snrRet.nii'];
    Vsnr = rmfield(Vsnr, 'pinfo');
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
            %Wp = [0.01/(fs/2) LPF/(fs/2)];
            %Ws = [0.8*Wp(1) 1.2*Wp(2)];
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
    dspTc2 = Y_s(:, roi_indx)';
    %% Pre-processed signal after LPF
    %%{
    Data_mat = zeros(im_size(1:3));
    for i = 1:nimgs
        Data_mat(keeps) = Y_s(i, :);
        Yimgs(:, :, :, i) = Data_mat;
    end
    Data_mat = zeros(im_size(1:3));
    snr_lpf = post_processMap(mean(Y_s)./std(Y_s));
    Data_mat(keeps) = snrMask.*snr_lpf;
    Vsnr = spm_vol(imgs{1});
    Vsnr.fname = [imgpath_roi filesep 'snrLpf.nii'];
    Vsnr = rmfield(Vsnr, 'pinfo');
    spm_write_vol(Vsnr, Data_mat);
    Xrrr = squeeze(Data_mat(:, :, dispSlice));
    disp('Writing out pre-processed images');
    if stat_4DFile == 0
        [p f e] = fileparts(imgs{1});
        ofuncs = spm_select('FPList', [p '/'], ['^' brick_img '.*\.nii$']);
        matlabbatch{1}.spm.util.cat.name = [p '/' brick_img '4D.nii'];
        matlabbatch{1}.spm.util.cat.vols = cellstr(ofuncs);
        matlabbatch{1}.spm.util.cat.dtype = 4;
        spm_jobman('run', matlabbatch);
        zpFile = matlabbatch{1}.spm.util.cat.name;
        delete([zpFile(1:end - 4) '.mat']);
        for c = 1:nimgs
            [p f e] = fileparts(imgs{c});
            V = spm_vol(imgs{c});
            spm_write_vol(V, squeeze(Yimgs(:, :, :, c)));
        end
        ofuncs = spm_select('FPList', [p '/'], ['^' brick_img '.*\.nii$']);
        matlabbatch{1}.spm.util.cat.name = [p '/prep_' brick_img '4D.nii'];
        matlabbatch{1}.spm.util.cat.vols = cellstr(ofuncs);
        matlabbatch{1}.spm.util.cat.dtype = 4;
        spm_jobman('run', matlabbatch);
        zpFile = matlabbatch{1}.spm.util.cat.name;
        delete([zpFile(1:end - 4) '.mat']);
    else
        for c = 1:nimgs
            [p f e] = fileparts(imgs{c});
            V = spm_vol(imgs{c});
            spm_write_vol(V, squeeze(Yimgs(:, :, :, c)));
        end
        %% Writting zipped 4D nifti file
        disp('Writting Preprocessed 4D Image File ...');
        clear matlabbatch
        ofuncs = spm_select('FPList', [p '/'], '^rest.*\.nii$');
        matlabbatch{1}.spm.util.cat.name = [p '/prep_rimg4D.nii'];
        matlabbatch{1}.spm.util.cat.vols = cellstr(ofuncs);
        matlabbatch{1}.spm.util.cat.dtype = 4;
        spm_jobman('run', matlabbatch);
        zpFile = matlabbatch{1}.spm.util.cat.name;
        delete([zpFile(1:end - 4) '.mat']);
        %delete([zpFile '.gz']);
    end
    P = squeeze(median(permute(Yimgs, [4 1 2 3])));
    V.fname = [p '/prep_medImg.nii'];
    spm_write_vol(V, P);
    Y_s = post_processMap(Y_s);
    %% Write filtered snr
    Data_mat = zeros(im_size(1:3));
    snr_preprocessed = mean(Y_s)./std(Y_s);
    snr_preprocessed = post_processMap(snr_preprocessed);
    Data_mat(keeps) = snrMask.*snr_preprocessed;
    Vsnr = spm_vol(imgs{1});
    Vsnr.fname = [imgpath_roi filesep 'snrPreprocessed.nii'];
    Vsnr = rmfield(Vsnr, 'pinfo');
    spm_write_vol(Vsnr, Data_mat);

    clear zmot zglob zhpf conf prin_c M VM coefs data Y1;
    B = zeros(nimgs, max(size(keeps)));
    Z = zeros(nimgs, max(size(keeps)));
    Data_mat = zeros(im_size(1:3));

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
    %% Signal power vs. frequency
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
    
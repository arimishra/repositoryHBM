%% spm_cnvt_viz.m
%
% Steady-state functional connectivity analysis toolbox for SPM5.  This
% function does vizualization of the connectivity map.
%
% Authors:
%     Arabinda Mishra
%     Baxter P. Rogers
%
% $Id: spm_cnvt.viz.m 1 2009-09-02 12:40:03Z mishra $

function spm_smclust(varargin)
    warning off;
    %% varargin arguments are
    %   1  roi_msk                     (ROI)
    %   2  n_grps                (Number of groups)
    %   3  spm_file              (Connectivity maps)
    %   4  net.trainParam.epochs (Connectivity threshold)
    %   5  somsv                 (Save clustered connectivity map flag)
    %   6  maptype               (Type of data handled, i.e., Connectivity map or time series data)
    %   7 extMask                (External Mask)

    maptype = varargin{1};
    spm_file = varargin{2};
    roi_msk = varargin{3};
    n_grps = varargin{5};
    extMask = varargin{7};
    [p1 f1 e1] = fileparts(roi_msk{1});
    
    Vout = spm_vol(roi_msk{1});
    V = spm_read_vols(Vout);
    im_size = size(V);
    ind = find(V > 0);
    [p f e] = fileparts(spm_file{1});
    Vout = spm_vol([p filesep f e(1:4)]);
    M = spm_read_vols(Vout);
    nimgs = size(M, 4);
    %% Self organizing map (SOM) parameterization
    if ~isempty(extMask{1})
        Vout = spm_vol(extMask{1});
        V = (spm_read_vols(Vout) > 0);
        keeps = find(squeeze(M(:, :, :, 1)).*V > 0);
    else
        V = ones(size(V));
        keeps = find(squeeze(M(:, :, :, 1)) > 0);
    end
    if maptype == 1
        Y_s = zeros(nimgs, length(keeps));
        for i = 1:nimgs
            X = squeeze(M(:, :, :, i));
            Y_s(i, :) = X(keeps);
        end
    else
        Y_s = zeros(length(ind), nimgs);
        for i = 1:nimgs
            X = squeeze(M(:, :, :, i));
            Y_s(:, i) = X(ind);
        end
    end
    
    Y_s = post_processMap(Y_s);
    Y_sM = Y_s;
    %% Reduce dimension of Connectivity Map include following lines
    %%{
    scaleFact = 3;
    if scaleFact > 1
        k = 1;
        for i = 1:scaleFact:length(keeps) - scaleFact
            Y_ss(:, k) = Y_s(:, i);
            k = k + 1;
        end
        Y_s = Y_ss;
    end
    %}
    net = selforgmap([1 n_grps]);
    net.trainParam.epochs = varargin{4};
    net = train(net, Y_s');
    a = net(Y_s');
    ac = vec2ind(a);
    data = zeros(1, im_size(1)*im_size(2)*im_size(3));
    data(ind) = ac;
    data = reshape(data, [im_size(1) im_size(2) im_size(3)]);
    con_map = Y_s;
    net_wt = net.iw{1, 1};

    %% Write clustered roi region map and SOM parameters
    Vout = spm_vol(roi_msk{1});
    if maptype == 1
        Vout.fname = [p1 filesep 'Fnc_Clst(' f1 num2str(n_grps) ').nii'];
    else
        Vout.fname = [p1 filesep 'FClst_Tseries(' f1 num2str(n_grps) ').nii'];
    end
    spm_write_vol(Vout, data);
    %% Include two lines below to save parameter files
    %% which are normally bulky
    %Vout.fname = [mpath filesep 'Fnc_Clst(' f num2str(n_grps) ')param.mat'];
    %save(Vout.fname, 'con_map', 'net_wt');
    somsv = varargin{5};

    %% Write the clustered connectivity maps and
    %% corresponding mean connectivity maps
    %{
    if maptype == 1
        kp = zeros(n_grps, length(ind));
        k = zeros(1, n_grps);
        for i = 1:max(size(ind))
            for j = 1:n_grps
                if data(ind(i)) == j
                    k(j) = k(j) + 1;
                    kp(j, k(j)) = i;
                end
            end
        end
        for i = 1:n_grps
            data = zeros(im_size);
            for j = 1:k(i)
                c = kp(i, j);
                X = zeros(1, im_size(1)*im_size(2)*im_size(3));
                X(keeps) = Y_sM(c, :);
                X = reshape(X, im_size(1:3));
                if somsv == 1
                    Vout.fname = [p1 filesep f(1:end - n_digit) num2str(n_grps) '_' num2str(i) ')_' f(end - n_digit + 1:end) '.nii'];
                    spm_write_vol(Vout, X);
                end
                data = data + X;
            end
            data = data/j;
            Vout.fname = [p1 filesep f(1) 'AvgCnvt' f(1) f1 'Seg' num2str(n_grps) '_' num2str(i) '.nii'];
            spm_write_vol(Vout, data);
            %data(keeps) = net_wt(i, :);
            %Vout.fname = [p1 filesep f(1:end - n_digit) '_Neuron' num2str(n_grps) '_' num2str(i) '.nii'];
            %spm_write_vol(Vout, data);
        end
    else
        data = dataM;
        X = squeeze(data(:, :, :, 10));
        imgthr = spm_antimode(X);
        keeps = find(X > imgthr/2.5);
        Y_ss = zeros(nimgs, length(keeps));
        for i = 1:nimgs
            X = squeeze(data(:, :, :, i));
            Y_s1(i, :) = X(keeps);
        end
        Y_s1 = post_processMap(Y_s1);
        %{
        for i = 1:n_grps
            keep = find(ac == i);
            R = mean(Y_s(keep, :))';
            X = zeros(1, im_size(1)*im_size(2)*im_size(3));
            R = R - mean(R);
            X(keeps) = atanh(corrcoef(Y_s1, R)) * sqrt(nimgs - 3);
            data = reshape(X, [im_size(1) im_size(2) im_size(3)]);
            %% Average connectivity maps
            Vout.fname = [p1 filesep 'Cnvty(Ts' f1 num2str(n_grps) '_' num2str(i) ')Avg.nii'];
            spm_write_vol(Vout, data);
            
            R = net_wt(i, :) - mean(net_wt(i, :));
            X(keeps) = atanh(corr(Y_s1, R')) * sqrt(nimgs - 3);
            data = reshape(X, [im_size(1) im_size(2) im_size(3)]);
            %% Neuron Connectivity maps
            Vout.fname = [p1 filesep 'CnvtyNrnTs(' f1 num2str(n_grps) '_' num2str(i) ').nii'];
            spm_write_vol(Vout, data);
        end
        %}
    end
    %}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
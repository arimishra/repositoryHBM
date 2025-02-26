%% spm_Dyn_cnvtOptmz.m
%
% Steady-state connectivity dynamics toolbox for SPM8/12.  This
% function does all the actual work.
% Author:
% Arabinda Mishra
% $Id: spm_stdcnvty.m 13 20018-11-27 10:40:03Z mishra $
% We can generate a visual report for each analysis, including slices of
% the connectivity map overlaid on the first functional image, and a
% histogram of the connectivity values.

    function spm_Dyn_cnvtOptmz(varargin)
        %%varargin arguments are
        %   1  imgs    (Image list)
        %   2  roi     (ROI file)

        imgs           = varargin{1};
        roi            = varargin{2};
        stat_4DFile = 0;
        [p f e] = fileparts(imgs{1});
        if e(1:3) == '.gz';
            gunzip([p '/' f e]);
            imgs{1} = [p '/' f];
            [p f e] = fileparts(imgs{1});
        end
        Vout = spm_vol(imgs{1});
        M = spm_read_vols(Vout);
        if sum(size(M) ~= 0) == 4
            gzip(imgs{1});
            [imgs] = fCnvty4_3Dconv(imgs{1});
            stat_4DFile = 1;
            zpFile = f;
            delete([p '/' f e]);
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
        %% Image input
        if ~isempty(imgs{1}) == 1
            nimgs = length(imgs);
            [imgpath img1file img1ext] = fileparts(imgs{1});
            if ~isempty(roi)
                if ~isempty(roi{1})
                    [imgpath_roi img1file_roi img1ext_roi] = fileparts(roi{1});
                end
            end
            Vimgs = spm_vol(strvcat(imgs));
            Yimgs = spm_read_vols(Vimgs);
            %% Image masking - use external mask if supplied; if not, calculate a
            %% histogram-based threshold
            %% imgthr should be imgthr for monkey data
            Vmask = spm_vol(roi{1});
            Ymask = spm_read_vols(Vmask);
            keeps = find(round(Ymask(:)) > 0);
            im_size = size(Ymask);

            %% Include line below if you want to see maps with all frames (without excluding the beginning volumes) 
            Y_s = zeros(nimgs, length(keeps));
            Data_mat = zeros(im_size(1:3));
            for i = 1:nimgs
                Data_mat = squeeze(Yimgs(:, :, :, i));
                Y_s(i, :) = Data_mat(keeps);
            end
            Y_s = post_processMap(Y_s);
            roi_indx = round(length(keeps)/3);
            
            slWin_width = 60;
            imgpath = imgpath_roi;
            stPnt = round(nimgs*rand) + 1;
            no_peak = 4;
            while (stPnt > nimgs - slWin_width),
                stPnt = stPnt - slWin_width;
            end
            while (stPnt < slWin_width),
                stPnt = stPnt + slWin_width;
            end
            r_del = 1;
            j = 1;
            YsT = Y_s(stPnt:stPnt + slWin_width - 1, :);
            while (r_del > 10^ - 10),
                for c = 1:nimgs - slWin_width
                    YsI = Y_s(c:c + slWin_width - 1, :);
                    muT = sum(YsT(:))/(length(keeps)*slWin_width);
                    muI = sum(YsI(:))/(length(keeps)*slWin_width);
                    stdP = sqrt(sum(sum((YsT - muT).^2)))*sqrt(sum(sum((YsI - muI).^2)));
                    r_val(c) = sum(sum((YsI - muI).*(YsT - muT)))/stdP;
                end

                r_ival(j, :) = r_val;
                for s = 1:no_peak
                    r_max = - 1;
                    rm_index = 0;
                    for i = 1:c
                        if r_val(i) > r_max
                            r_max = r_val(i);
                            rm_index = i;
                        end
                    end
                    r_peak(j, s) = rm_index;
                    r_val(rm_index) = - 1;
                end
                YsT = zeros(size(YsT));
                for s = 1:no_peak
                    YsT = YsT + Y_s(r_peak(j, s):r_peak(j, s) + slWin_width - 1, :);
                end
                YsT = YsT/no_peak;
                
                if j > 1
                    X = r_ival(j, :);
                    Y = r_ival(j - 1, :);
                    r_del = sum(X(r_peak(j, :)) - Y(r_peak(j - 1, :)));
                end
                j = j + 1;
            end
            
            Data_mat = zeros(im_size(1:3));
            for c = 1:nimgs - slWin_width
                r_vl = corr(Y_s(c:c + slWin_width - 1, :), YsT);
                for k = 1:length(keeps)
                    rVal(k) = r_vl(k, k);
                end
                Data_mat(keeps) = rVal;
                r_t(c) = sum(Data_mat(:))/length(keeps);
                %{
                figure;
                imagesc(squeeze(Data_mat(:, :, 3)));
                axis off;
                colormap(jet);
                colorbar;
                %}
                Vout = spm_vol(imgs{1});
                Vout.fname = [imgpath filesep 'TempDyn_' num2str(c, '%03.0f') '.nii'];
                spm_write_vol(Vout, Data_mat);
            end
            %figure;
            %plot(1.5*(1:50), r_t(1:50))
            %}
        end
        %% Delete 3D volumes
        %%{
        brick_img = 'TempDyn';
        [p f e] = fileparts(imgs{1});
        ofuncs = spm_select('FPList', [p '/'], ['^' brick_img '.*\.nii$']);
        matlabbatch{1}.spm.util.cat.name = [p '/' brick_img '4D.nii'];
        matlabbatch{1}.spm.util.cat.vols = cellstr(ofuncs);
        matlabbatch{1}.spm.util.cat.dtype = 4;
        spm_jobman('run', matlabbatch);
        gzip(matlabbatch{1}.spm.util.cat.name);
        zpFile = matlabbatch{1}.spm.util.cat.name;
        delete([zpFile(1:end - 4) '.mat']);
        delete(zpFile);
        if stat_4DFile == 1
            [p f e] = fileparts(imgs{1});
            for c = 1:nimgs
                delete([p '/rest' num2str(c, '%04.0f') '.nii']);
            end
        end
        for c = 1:nimgs - slWin_width
            delete([p '/TempDyn_' num2str(c, '%03.0f') '.nii']);
        end
        %}
        %{
        figure;
        %% Display the time course of four representative voxels.
        subplot(2, 2, 1), plot(dspTc, 'r'), hold on;
        xlabel('Ex. Time Course: RAW DATA');
        title('REPRESENTATIVE VOXEL TIME COURSE');

        subplot(2, 2, 2), plot(dspTc1, 'g'), hold on;
        xlabel('Ex. Time Course: REGRESSED DATA');

        subplot(2, 2, 3), plot(dspTc2, 'b'), hold on;
        xlabel('Ex. Time Course: LOW PASS FILTERED');

        subplot(2, 2, 4), plot(dspTc, 'r'), hold on;
        subplot(2, 2, 4), plot(dspTc1, 'g'), hold on;
        subplot(2, 2, 4), plot(dspTc2, 'b'), hold on;
        xlabel('TIME COURSE OVERLAID');
        close all;
        %}
        
%% Convert fdf to Nifti
%% Program revised by Arabinda Mishra June 10, 2019
%% 1st step of pre-processing the raw Varian input data in .fdf or available .nii 4D data
%% Run time close to 2min on 64bit mac of 300 volumes
clear all;
close all;
warning off;
    setenv('FSLDIR', '/usr/local/fsl');
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
    disp('Select any .FDF/.nii or zipped nifti file to Process ....');
    [filename, outputPath] = uigetfile('*.*', 'Select image file/files');
    imType = 'Fnc';
    %% Number of volumes want to delete off the begining
    delVol = 0;
    if filename(end - 2:end) == 'fdf'
        c = 1;
        while (exist([outputPath filename(1:5) num2str(c, '%03.0f') filename(9:end)]))
            c = c - 1;
        end
        c = c + 1;
        while (exist([outputPath filename(1:5) num2str(c, '%03.0f') filename(9:end)]))
            c = c + 1;
        end
        num_slice = c - 1;
        c = str2num(filename(14:16));
        while (exist([outputPath filename(1:13) num2str(c, '%03.0f') filename(17:end)]))
            c = c - 1;
        end
        c = c + 1;
        while (exist([outputPath filename(1:13) num2str(c, '%03.0f') filename(17:end)]))
            c = c + 1;
        end
        c = c - 1;
        num_vol = c;
        if num_vol > 1
            for c = delVol + 1:num_vol
                for i = 1:num_slice
                    slice_name = [outputPath filename(1:5) num2str(i, '%03.0f') filename(9:13) num2str(c, '%03.0f') filename(17:end)];
                    slice = vuOpenImage_nii(slice_name);
                    im.Data(:, :, i) = slice.Data;
                end
                im.Spc = slice.Spc;
                im.Origin = slice.Origin;
                im.Dims = [slice.Dims num_slice];
                im.Orientation = slice.Orientation;
                write_file_name = [outputPath 'image', num2str(c, '%03.0f') '.nii'];
                vuWriteImage_nii(im, write_file_name); 
            end
        elseif num_vol == 1
            for i = 1:num_slice
                slice_name = [outputPath filename(1:5) num2str(i, '%03.0f') filename(9:13) num2str(c, '%03.0f') filename(17:end)];
                slice = vuOpenImage_nii(slice_name);
                im.Data(:, :, i) = slice.Data;
            end
            im.Spc = slice.Spc;
            im.Origin = slice.Origin;
            im.Dims = [slice.Dims num_slice];
            im.Orientation = slice.Orientation;
            write_file_name = [outputPath 'anatImage.nii'];
            vuWriteImage_nii(im, write_file_name);
            disp(['Writting Anatomic Image into the Folder' outputPath]);
            imType = 'Ant';
        end
    elseif filename(end - 2:end) == 'nii'
        img4_3Dconv([outputPath filename])
    elseif filename(end - 1:end) == 'gz'
        img4_3Dconv([outputPath filename])
    else
        disp('CHECK IMAGE FORMAT ....')
    end
    if imType == 'Fnc'
        ofuncs = spm_select('FPList', outputPath, '^image.*\.nii$');
        spm_realign(ofuncs);
        if filename(end - 2:end) == 'fdf'
            rp_mot = load([outputPath 'rp_image' num2str(delVol + 1, '%03.0f') '.txt']);
            copyfile([outputPath 'rp_image' num2str(delVol + 1, '%03.0f') '.txt'], [outputPath 'rp_mot3D.txt']);
            delete([outputPath 'rp_image' num2str(delVol + 1, '%03.0f') '.txt']);
        else
            rp_mot = load([outputPath 'rp_image001.txt']);
            copyfile([outputPath 'rp_image001.txt'], [outputPath 'rp_mot3D.txt']);
            delete([outputPath 'rp_image001.txt']);
        end
        rp_mot(:, 4:6) = 180*rp_mot(:, 4:6)/pi;
        ofuncs = spm_select('FPList', outputPath, '^image.*\.nii$');
        matlabbatch{1}.spm.util.cat.name = [outputPath 'img4D.nii'];
        matlabbatch{1}.spm.util.cat.vols = cellstr(ofuncs);
        matlabbatch{1}.spm.util.cat.dtype = 4;
        spm_jobman('run', matlabbatch);
        zpF = matlabbatch{1}.spm.util.cat.name;
        delete([zpF(1:end - 4) '.mat']);
        if max(abs(max(rp_mot))) > 2
            disp([subName{p} ' : Run -' num2str(c) ' may have motion artifacts...']);
            figure;
            subplot(2, 1, 1), plot(rp_mot(:, 4), 'r', 'LineWidth', 1); hold on;
            subplot(2, 1, 1), plot(rp_mot(:, 5), 'b', 'LineWidth', 1); hold on;
            subplot(2, 1, 1), plot(rp_mot(:, 6), 'g', 'LineWidth', 1);
            title(['Motion Parameters (WARNING.... too much motion .. ' subName{p} ' -Run: ' num2str(c)]);
            xlabel('Number of dynamics');
            ylabel('Rotation (degrees)');
            subplot(2, 1, 2), plot(rp_mot(:, 1), 'r', 'LineWidth', 1); hold on;
            subplot(2, 1, 2), plot(rp_mot(:, 2), 'b', 'LineWidth', 1); hold on;
            subplot(2, 1, 2), plot(rp_mot(:, 3), 'g');
            xlabel('Number of dynamics');
            ylabel('Displacement (mm)');
        else
            figure;
            subplot(2, 1, 1), plot(rp_mot(:, 4), 'r', 'LineWidth', 1); hold on;
            subplot(2, 1, 1), plot(rp_mot(:, 5), 'b', 'LineWidth', 1); hold on;
            subplot(2, 1, 1), plot(rp_mot(:, 6), 'g', 'LineWidth', 1);
            title('Motion Parameters');
            xlabel('Number of dynamics');
            ylabel('Rotation(digrees)');
            subplot(2, 1, 2), plot(rp_mot(:, 1), 'r', 'LineWidth', 1); hold on;
            subplot(2, 1, 2), plot(rp_mot(:, 2), 'b', 'LineWidth', 1); hold on;
            subplot(2, 1, 2), plot(rp_mot(:, 3), 'g');
            xlabel('Number of dynamics');
            ylabel('Displacement(mm)');
        end
        im2b_SLTcorrect = [outputPath 'img4D.nii'];
        SLTCorrectImage = [outputPath 'rimg4D.nii'];
        disp('In process ...');
        %% down or odd (--down, --odd) for descending and interleaved acquisition for SLICE TIME CORRECTION ...
        %% Default TR = 3s
        inputFSLparms = ['/usr/local/fsl/bin/slicetimer -i ' im2b_SLTcorrect ' -r 3 --odd -o ' SLTCorrectImage];
        unix(inputFSLparms);
        gunzip([outputPath 'rimg4D.nii.gz']);
        delete([outputPath 'rimg4D.nii.gz']);
        Vout = spm_vol(SLTCorrectImage);
        Yimgs = spm_read_vols(Vout);
        M = squeeze(median(permute(Yimgs, [4 1 2 3])));
        Vout = spm_vol([outputPath 'image001.nii']);
        Vout.fname = [outputPath 'med_img.nii'];
        spm_write_vol(Vout, M);
        for c = 1:length(ofuncs) + 100
            delete([outputPath 'image' num2str(c, '%03.0f') '.nii']);
        end
        disp(['Writting EPI image into Folder' outputPath]);
    end
    
    
    
    

   


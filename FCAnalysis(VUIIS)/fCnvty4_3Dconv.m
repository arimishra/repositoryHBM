    function [imgs] = fCnvty4_3Dconv(file_Name)
        Vout = spm_vol(file_Name);
        M = spm_read_vols(Vout);
        [pathname f e] = fileparts(file_Name);
        for i = 1:size(M, 4)
            if i == 1
                V.fname = [pathname '/rest' num2str(i, '%04.0f') '.nii'];
                V.mat = Vout(i).mat;
                V.dim = Vout(i).dim;
                V.dt = Vout(i).dt;
                V.pinfo = Vout(i).pinfo;
                V.n = Vout(i).n;
                V.descrip = '3D data';
                V.private.dat.fname = V.fname;
                V.private.dat.dim = Vout(i).dim;
                V.private.dat.dtype = Vout(i).private.dat.dtype;
                V.private.dat.offset = Vout(i).private.dat.offset;
                V.private.dat.scl_slope = Vout(i).private.dat.scl_slope;
                V.private.dat.scl_inter = Vout(i).private.dat.scl_inter;
                V.private.mat = Vout(i).private.mat;
                V.private.mat_intent = Vout(i).private.mat_intent;
                V.private.mat0 = Vout(i).private.mat0;
                V.private.mat0_intent = Vout(i).private.mat0_intent;
                V.private.descrip = '3D data';
                spm_write_vol(V, squeeze(M(:, :, :, i)));
            else
                V = spm_vol([pathname '/rest' num2str(1, '%04.0f') '.nii']);
                V.fname = [pathname '/rest' num2str(i, '%04.0f') '.nii'];
                spm_write_vol(V, squeeze(M(:, :, :, i)));
            end
        end
        for c = 1:size(M, 4)
            Vout = spm_vol([pathname '/rest' num2str(c, '%04.0f') '.nii']);
            M = spm_read_vols(Vout);
            Vout = spm_vol([pathname '/rest0001.nii']);
            Vout.fname = [pathname '/rest' num2str(c, '%04.0f') '.nii'];
            spm_write_vol(Vout, M);
            imgs{c} = [pathname '/rest' num2str(c, '%04.0f') '.nii'];
        end
        %% Special case for un interpolated rat images
        %{
        if size(M, 1) == 96
            [x y X Y] = interpAnatFnc(zeros(96), zeros(128));
            mF = 128/96;
            num_slice = size(M, 3);
            for c = 1:size(M, 4)
                Mp = zeros([128 128 num_slice]);
                Vout = spm_vol([pathname '/rest' num2str(c, '%04.0f') '.nii']);
                Vp = spm_read_vols(Vout);
                for i = 1:num_slice
                    Mp(:, :, i) = interp2(x, y, squeeze(Vp(:, :, i)), X, Y, 'linear');
                end
                Vout.dim = [128 128 num_slice];
                Vout.mat(1, 1) = Vout.mat(1, 1)/mF;
                Vout.mat(2, 2) = Vout.mat(2, 2)/mF;
                Vout.private.dat.dim = Vout.dim;
                Vout.private.mat = Vout.mat;
                spm_write_vol(Vout, Mp);
            end
        end
        %}
        
        
        
        
        
        
        
        
%% 4D23D converter
function img4_3Dconv(fileName)
        if fileName(1:end - 2:end) == 'gz'
            gunzip(fileName);
            Vout = spm_vol(fileName(1:end - 3));
        else
            Vout = spm_vol(fileName);
        end
        M = spm_read_vols(Vout);
        [pathname f e] = fileparts(fileName);
        for i = 1:size(M, 4)
            if i == 1
                V.fname = [pathname '/image' num2str(i, '%03.0f') '.nii'];
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
                V = spm_vol([pathname '/image' num2str(1, '%03.0f') '.nii']);
                V.fname = [pathname '/image' num2str(i, '%03.0f') '.nii'];
                spm_write_vol(V, squeeze(M(:, :, :, i)));
            end
        end
        if fileName(1:end - 2:end) == 'gz'
            delete(fileName(1:end - 3));
        end
    end
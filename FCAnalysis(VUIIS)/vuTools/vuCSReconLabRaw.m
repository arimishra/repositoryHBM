function I = vuCSReconLabRaw (filename, opts)
%%VUCSRECONLABRAW Reconstruct CS acquisitions on Philips scanners from LAB/RAW data.

[data info] = loadLabRaw(filename);
data = squeeze(data);
[N1 N2 N3] = size(data);

% ky values
ky = double( info.labels.RTopOffset.vals(info.labels_row_index_array(:)) );
idx = find(ky > 2^15 - 1);
ky(idx) = ky(idx) - 2^16;

% kz values
kz = double( info.labels.RRInterval.vals(info.labels_row_index_array(:)) );
idx = find(kz > 2^15 - 1);
kz(idx) = kz(idx) - 2^16;

% Actual data dimensions
Nkx = N1/2;
Nky = 2*max(abs(ky(:)));
Nkz = 2*max(abs(kz(:)));

% Place readouts into proper ky, kz locations
data_sorted = zeros(2*Nkx,Nky,Nkz);
M = false(2*Nkx,Nky,Nkz);
for idx2 = 1:N2
	for idx3 = 1:N3
		idxlabel = info.labels_row_index_array(idx2,idx3);
		
		thisky = double( info.labels.RTopOffset.vals(idxlabel) );
		if thisky > 2^15 - 1
			thisky = thisky - 2^16;
		end
		
		thiskz = double( info.labels.RRInterval.vals(idxlabel) );
		if thiskz > 2^15 - 1
			thiskz = thiskz - 2^16;
		end
		
		data_sorted(:, thisky+Nky/2+1, thiskz+Nkz/2+1) = data(:,idx2,idx3);
		M(:, thisky+Nky/2+1, thiskz+Nkz/2+1) = true;
	end
end

% CS Recon
I = vuCSReconCart(fftshift(data_sorted),fftshift(M),opts);
I = fftshift(I,1);

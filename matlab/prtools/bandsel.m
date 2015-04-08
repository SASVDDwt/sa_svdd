%BANDSEL Selection of bands from object images
%
%   B = BANDSEL(A,J)
%   W = BANDSEL([],J)
%   B = A*BANDSEL([],J)
%
% INPUT
%   A    Dataset or datafile with multi-band object images
%   J    Bands to be selected
%
% OUTPUT
%   W    Mapping performing the band selection
%   B    Dataset with selected bands (oredered according to J)
%
% DESCRIPTION
% If the objects in a dataset or datafile A are multi-band images, e.g. RGB
% images, or the result of IM_PATCH, then the featsize of A is [M,N,L] for
% for L bands of an M x N images. This routine makes a selection J out of
% L. The routine BAND2OBJ may be used BANDS vertically as separate objects.
%
% SEE ALSO
% DATASETS, DATAFILES, IM2OBJ, IM_PATCH, BAND2OBJ

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = bandsel(a,J)

mapname = 'Band Selection';
if nargin < 2, J = 1; end
if nargin < 1 | isempty(a)
	w = mapping(mfilename,'fixed',{J});
	w = setname(w,mapname);
elseif isdatafile(a)
	isobjim(a);
	fsize = getfeatsize(a);
	v = mapping(mfilename,'fixed',{[]},[],fsize,[fsize(1:2) length(J)]);
	% make the mappings identical, i.e. independent of v.data
	% and store the data in the ident field
	% This will enable vertical concatenation of datafiles
	a = setident(a,J,'bandsel');
	v = setdata(v,[]);
	w = addpostproc(a,v);
elseif isdataset(a)
	m = size(a,1);
	if isempty(J)  % J is stored in a.ident
		J = getident(a,'bandsel');
	else
		J = repmat(J(:)',m,1);
	end
	isobjim(a);
	fsize = getfeatsize(a);
	if length(fsize) < 3
		error('No image bands defined for dataset')
	end
	if any(J(:) > fsize(3))
		error('Wrong index for image bands')
	end
	k = prod(fsize(1:2));                           % size of a band
	L = repmat((J(:,1)-1)*k,1,k)+repmat([1:k],m,1); % indices of band_1 per object
	if size(J,2) > 1                                % concatenate in case of multiple band selection
		for j=2:size(J,2)
			LL = repmat((J(:,j)-1)*k,1,k)+repmat([1:k],m,1); % indices of band_j per object
			L = [L LL];                                % indices for all bands to be selected
		end
	end
	adata = getdata(a);                            % Let us do the selection on the data
	bdata = zeros(m,k*size(J,2));
	for i=1:m                                      % can this be done faster?
		bdata(i,:) = adata(i,L(i,:));
	end
	b = setdat(a,bdata);                           % store the data in a dataset with all information of a
	w = setfeatsize(b,[fsize(1:2) size(J,2)]); 
else
	error('Illegal command')
end
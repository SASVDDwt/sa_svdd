%BAND2OBJ Mapping image bands to objects
%
%		B = BAND2OBJ(A,N)
%   W = BAND2OBJ([],N)
%   B = A*W
%
% INPUT
%   A   Dataset or datafile with multiband image objects.
%   N   Number of successive bands to be combined in an object.
%       The number of image bands in A should be multiple of N.
%       Default N = 1.
%
% OUTPUT
%   W    Mapping performing the band selection
%   B   Output dataset or datafile.
%
% DESCRIPTION
% If the objects in a dataset or datafile A are multi-band images, e.g. RGB
% images, or the result of IM_PATCH, then the featsize of A is [C,R,L] for
% for L bands of an C x R images. This routine combines sets of N successive 
% bands as separate objects. The total number of objects is thereby
% enlarged by a factor L/N. All information of the constituting objects
% like labels, is copied to the newly created objects.
%
% SEE ALSO
% DATASETS, DATAFILES, BANDSEL, IM2OBJ, DATA2IM

function b = band2obj(a,n)

mapname = 'Bands to Objects';
if nargin < 2, n = 1; end
if nargin < 1 | isempty(a)
	b = mapping(mfilename,'fixed',{n});
	b = setname(b,mapname);
else
	isvaldfile(a);
	isobjim(a);
	k = getfeatsize(a,3);
	if n*round(k/n) ~= k
		error('Number of images bands should be multiple of N')
	end
	b = bandsel(a,1:n);
	for j=2:round(k/n)
		b = [b; bandsel(a,(j-1)*n+1:j*n)];
	end
end

return
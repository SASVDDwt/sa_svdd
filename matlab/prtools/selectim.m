%SELECTIM Select one or more images in multiband image or dataset
%
%	B = SELECTIM(A,N)
% A = A*SELECTIM([],N)
%
% INPUT
%   A          Multiband image or dataset containing multiband images
%   N          Vector or scalar pointing to desired images
%
% OUTPUT
%   B          New, reduced, multiband image or dataset

function b = selectim(a,n)

	if nargin < 2, n = 1; end
	if nargin < 1 | isempty(a)
		b = mapping(mfilename,'fixed',{n})
		b = setname(b,'SelectImage');
		return
	end
	
	if isdataset(a)
		im = data2im(a);
		if ndims(im) ~= 3
			error('3D images expected')
		end
		fsize = size(im);
		if length(n) == 1
			fsize = fsize(1:2);
		else
			fsize(3) = length(n);
		end
		im = im(:,:,n);
		if isobjim(a)
			b = setdat(a,im2obj(im,fsize));
			b = setfeatsize(b,fsize);
		else
			b = im2feat(im);
		end
	elseif isdatafile(a)
		b = a*filtm([],'selectim',n);
	elseif isa(a,'double')
		if ndims(a) ~= 3
			error('3D images expected')
		end
		b = a(:,:,n);
	else
		error('Unexpected datatype')
	end
	
return
		

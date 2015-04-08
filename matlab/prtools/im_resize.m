%IM_RESIZE Mapping for resizing object images in datasets and datafiles
%
%  B = IM_RESIZE(A,SIZE,METHOD)
%  B = A*IM_RESIZE([],SIZE,METHOD)
%
% INPUT
%  A       Dataset or datafile
%  SIZE    Desired size
%  METHOD  Method, see IMRESIZE
%
% OUTPUT
%  B       Dataset or datafile
%
% DESCRIPTION
% The objects stored as images in the dataset or datafile A are resized
% using the IMRESIZE command. Default METHOD is 'nearest' (nearest neighbor
% interpolation. In SIZE the desired output size has to be stored. Note
% that for proper use in PRTools third size parameter of multi-band images,
% e.g. 3 for RGB images, has to be supplied. 
%
% SEE ALSO
% MAPPINGS, DATASETS, DATAFILES, IM2OBJ, DATA2IM 

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

%DXD 24-8-2007
%  I rewrote a part of this function. Now there are default values
%  given, a bug is removed, and the identation is correct again.
function b = im_resize(a,imsize,method)

if nargin < 3 | isempty(method)
	method = 'nearest';
end
if nargin<2 | isempty(imsize)
	imsize = [16 16];
end
if nargin<1
	a = [];
end

if isempty(a)
	b = mapping(mfilename,'fixed',{imsize,method});
% 	b = setname(b,'image resize');
% 	if length(imsize) == 1
% 		error('Mappings should have fixed output size')
% 	end
% 	b = setsize_out(b,imsize);
elseif isa(a,'dataset') % allows datafiles too
	isobjim(a);
	%b = filtim(a,mfilename,{imsize,method},imsize);
	b = filtim(a,mfilename,{imsize,method});
elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
	b = double(a);
	if length(imsize) > 1
		b = imresize(b,imsize(1:2),method);
	elseif imsize > 1 & round(imsize) == imsize
		b = imresize(a,[imsize imsize],method);
	else
		[m,n] = size(a);
		b = imresize(a,round(imsize*[m,n]),method);
	end
else
	class(a)
	error('Datatype not supported')
end



%FILTIM Mapping to filter multiband image objects in datasets and datafiles
%
%    B = FILTIM(A,FILTER_COMMAND,{PAR1,PAR2,....},SIZE)
%    B = A*FILTIM([],FILTER_COMMAND,{PAR1,PAR2,....},SIZE)
%
% INPUT
%    A               Dataset or datafile with image objects
%    FILTER_COMMAND  String with function name
%    {PAR1, ...  }   Cell array with optional parameters to FILTER_COMMAND
%    SIZE            Output size of the mapping (default: input size)
%
% OUTPUT
%    B               Dataset containing images processed by FILTER_COMMAND,
%                    image by image. Image bands are processes separately.
%
% DESCRIPTION
% For each object stored in A a filter operation is performed as
%
%    OBJECT_OUT = FILTER_COMMAND(OBJECT_IN,PAR1,PAR2,....)
%
% The results are collected and stored in B. In case A (and thereby B) is
% a datafile, execution is postponed until conversion into a dataset, or a
% call to SAVEDATAFILE.
%
% The difference between FILTIM and the similar command FILTM is that
% FILTIM is aware of the image structure and tests on it.
%
% EXAMPLE
% a = delft_images; b = a(120 121 131 230)*col2gray
% e = b*filtim([],'fft2')*filtim([],'abs')*filtim([],'fftshift');
% figure; show(e); figure; show((1+e)*filtim([],'log')); 
%
% SEE ALSO
% DATASETS, DATAFILES, IM2OBJ, DATA2IM, IM2FEAT, DATGAUSS, DATFILT, FILTM

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = filtim(a,command,pars,outsize)

	prtrace(mfilename);

	if nargin < 4, outsize = []; end
	if nargin < 3, pars = {}; end
	if nargin < 2
		error('No command given')
	end
	if ~iscell(pars), pars = {pars}; end
	
	mapname = 'dataset image processing';
	
	if isempty(a)                    % no data, so just mapping definition
		b = mapping(mfilename,'fixed',{command,pars});
		if ~isempty(outsize)
			b = setsize_out(b,outsize);
		end
		b = setname(b,mapname);
		return
	end

  if isdatafile(a)                 
		
		% for datafiles filters will be stored
		isobjim(a);      
				
    if isempty(getpostproc(a))     % as preprocessing (if no postproc defined)
      %b = addpreproc(a,mfilename,{command pars},outsize);
      b = addpreproc(a,command,pars,outsize);
    else                           % or as mapping as postprocessing
		  v = mapping(mfilename,'fixed',{command,pars});
		  if ~isempty(outsize)
			  v = setsize_out(v,outsize);
		  end
		  v = setname(v,mapname);
		  b = addpostproc(a,v);
		end
		if isempty(outsize)
			outsize = getfeatsize(b);
			b = setfeatsize(b,outsize);
		end
% 		if ~isempty(outsize)         % dont do this. It generates calls to
% 			b = setfeatsize(b,outsize);% readdatafile that may be bad for the
% 		end                          % cell datafile type
		return
    
	elseif isdataset(a)
   
		% convert to image and process
		isobjim(a);
    m = size(a,1);
    d = data2im(a);
		%out = feval(mfilename,d,command,pars{:}); % effectively jumps to double
		out = feval(mfilename,d,command,pars); % effectively jumps to double
		fsize = size(out);		
		if m > 1
			fsize = fsize(1:3);
			if fsize(3) == 1
				fsize = fsize(1:2);
			end
		end
    % store processed images in dataset
		b = setdata(a,im2obj(out,fsize));
		b = setfeatsize(b,fsize);
		
	else % double
		
		% make imsize 4D: horz*vert*bands*objects
    imsize = size(a);
		if length(imsize) == 2
			imsize = [imsize 1 1];
		elseif length(imsize) == 3
			imsize = [imsize 1];
		end
    
		% find size of first image
	  first = feval(command,getsubim(a,1,1),pars{:});
		if all(imsize(3:4) == 1) % preserve DipLib format for the time being
			b = first;
			return
		end
    first = double(first); % for Dip_Image users
		outsize = size(first);
		if length(outsize) == 3 % single subimage generates multiple bands !!
			b = repmat(first(:,:,1),[1 1 imsize(3)*outsize(3) imsize(4)]);
		else
			b = repmat(first,[1 1 imsize(3:4)]);
		end
    % process all other images
		nn = imsize(3)*imsize(4);
		s = sprintf('Filtering %i images',nn);
		prwaitbar(nn,s);
		for i = 1:imsize(3)
			for j = 1:imsize(4)
				prwaitbar(nn,(i-1)*imsize(4)+j);
				ima = double(feval(command,getsubim(a,i,j),pars{:}));
				if (any(outsize ~= size(ima)))  % check size
					error('All image sizes should be the same')
		  	end
				if length(outsize) == 2 % simple case: image_in --> image_out
					b(:,:,i,j) = ima;
				else                    % advanced: image_in --> bands out
					b(:,:,i:imsize(3):end,j) = ima;
				end
			end
		end
		prwaitbar(0);
  end
  
return

function asub = getsubim(a,i,j)

% needed as Dip_Image cannot handle 4d subsript if data is 2d or 3d

n = ndims(a);

if n == 2
	if i ~= 1 | j ~= 1
		error('Illegal image requested')
	end
	asub = a;
elseif n == 3
	if j ~= 1
		error('Illegal image requested')
	end
	asub = a(:,:,i);
elseif n == 4
	asub = a(:,:,i,j);
else
	error('Incorrect image supplied')
end
		

 

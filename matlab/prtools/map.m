%MAP Map a dataset, train a mapping or classifier, or combine mappings
%
%	B = MAP(A,W) or B = A*W
%
% Maps a dataset A by a fixed or trained mapping (or classifier) W, generating
% a new dataset B. This is done object by object. So B has as many objects
% (rows) as A. The number of features of B is determined by W. All dataset
% fields of A are copied to B, except the feature labels. These are defined
% by the labels stored in W.
%
%	V = MAP(A,W) or B = A*W
%
% If W is an untrained mapping (or classifier), it is trained by the dataset A.
% The resulting trained mapping (or classifier) is stored in V.
%
%	V = MAP(W1,W2) or V = W1*W2
%
% The two mappings W1 and W2 are combined sequentially. See SEQUENTIAL for
% a description. The resulting combination is stored in V.
%
% See also DATASETS, MAPPINGS, SEQUENTIAL

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: map.m,v 1.6 2007/06/26 10:30:01 duin Exp $

function [d, varargout] = map(a,b)
		prtrace(mfilename);
    
varargout = repmat({[]},[1, max((nargout-1),0)]);

[ma,ka] = size(a);
[mb,kb] = size(b);

if iscell(a) | iscell(b)
	if (iscell(a) & min([ma,ka]) ~= 1) | (iscell(b)& min([mb,kb]) ~= 1)
		error('Only one-dimensional cell arrays are supported')
	end

	if iscell(a) & ~iscell(b)
		d = cell(size(a));
		for j=1:length(a(:))
			d{j} = map(a{j},b);
		end
	elseif ~iscell(a) & iscell(b)
		d = cell(size(b));
		for j=1:length(b(:))
			d{j} = map(a,b{j});
		end
	else
		d = cell(length(a),length(b));
		for i = 1:length(a)
		for j = 1:length(b)
			d{i,j} = map(a{i},b{j});
		end
		end
	end
	return
end

if all([ma,mb,ka,kb] ~= 0) & ~isempty(a) & ~isempty(b) & ka ~= mb
	error('Output size of first argument should match input size of second')
end

if isa(a,'mapping') & isa(b,'mapping')
	
  if isempty(b) % empty mappings are treated as unity mappings
    d = a;
  elseif istrained(a) & isaffine(a) & istrained(b) & isaffine(b)
						% combine affine mappings
		d = affine(a,b);
	else
	  d = sequential(a,b);
	end
	
elseif isa(a,'dataset') | isa(a,'datafile') | isa(a,'double')
  
	if ~isa(b,'mapping')
		error('Second argument should be mapping or classifier')
	end

  if isempty(b) % treat empty mappings as unity mappings
    d = a;
    return
		% produce empty mapping (i.e. unity mapping)
		% if input data is empty
  elseif isempty(a)
		d = mapping([]);
    return
	  % handle scalar * mapping by .* (see times.m)
  elseif isa(a,'double') & ka == 1 & ma == 1
		d = a.*b;
		return; 
	end


	mapp = getmapping_file(b);

	if isuntrained(b)
		pars = +b;
		if iscombiner(feval(mapp)) % sequentiall, parallel and stacked need
                               % special treatment as untrained combiners
      
      % matlab 5 cannot handle the case [d, varargout{:}] = feval when varargout is empty
      % because of this we have next piece of code
			if isempty(varargout) 
        d = feval(mapp,a,b);   
      else    
        [d, varargout{:}] = feval(mapp,a,b); 
      end 
		else
			if isempty(varargout)
        d = feval(mapp,a,pars{:});
      else    
        [d, varargout{:}] = feval(mapp,a,pars{:});
      end 
		end
		if ~isa(d,'mapping')
			error('Training an untrained classifier should produce a mapping')
		end
		if getout_conv(b) > 1
			d = d*classc;
		end
		name = getname(b);
		if ~isempty(name)
			d = setname(d,name);
		end

	elseif isdatafile(a) & istrained(b)
		d = addpostproc(a,{b}); % just add mapping to postprocesing and execute later
		
	elseif isdatafile(a)
		try  % try whether this is a mapping that knows how to handle a datafile
			pars = getdata(b); % parameters supplied in fixed mapping definition
			if nargout > 0
				if isempty(varargout)
					d = feval(mapp,a,pars{:});
				else    
					[d, varargout{:}] = feval(mapp,a,pars{:});
     		end 
			else
				feval(mapp,a,pars{:})
				return
			end
		catch
      d = addpostproc(a,{b}); % just add mapping to postprocesing and execute later
      return
		end
  elseif isfixed(b) & isparallel(b)
    d = parallel(a,b);
	elseif isfixed(b) | iscombiner(b)
		pars = getdata(b); % parameters supplied in fixed mapping definition
		if nargout > 0
			if isempty(varargout)
				s_in = getsize_in(b);
				if prod(s_in) > 0
					a = setfeatsize(a,getsize_in(b)); % needed to set object images, sometimes
				end 
				d = feval(mapp,a,pars{:});  % Damn! sequential mappings are split! Solve there!
      else    
        [d, varargout{:}] = feval(mapp,a,pars{:});
      end 
		else
			feval(mapp,a,pars{:})
			return
		end

	elseif istrained(b)
		if ~isdataset(a) 
			a = dataset(a);
		end
		if isempty(varargout)
			fsize = getsize_in(b);
			if any(fsize~=0)
				a = setfeatsize(a,getsize_in(b)); % needed to set object images, sometimes
			end
			d = feval(mapp,a,b);
    else    
      [d, varargout{:}] = feval(mapp,a,b);
    end 
		if ~isreal(+d)
			prwarning(2,'Complex values appeared in dataset');
		end
		if isdataset(d)
			d = setcost(d,b.cost);
			% see if we have reasonable data in the dataset
	
		end

	else
		error(['Unknown mapping type: ' getmapping_type(b)])
	end
	
  if isdataset(d) 
	
			% we assume that just a basic dataset is returned, 
			% including setting of feature labels, but that scaling
			% and outputconversion still have to be done.

			% scaling
		v = getscale(b);
		if length(v) > 1, v = repmat(v(:)',ma,1); end
		d = v.*d;
			% outputconversion
		switch 	getout_conv(b);
		case 1  % SIGM output
			if size(d,2) == 1
				d = [d -d]; % obviously still single output discriminant
				d = setfeatlab(d,d.featlab(1:2,:));
			end             
			d = sigm(d);
		case 2  % NORMM output
			if size(d,2) == 1
				d = [d 1-d]; % obviously still single output discriminant
				d = setfeatlab(d,d.featlab(1:2,:));
			end             % needs conversion to two-classes before normm
			d = normm(d);
		case 3  % SIGM and NORMM output
			if size(d,2) == 1
				d = [d -d]; % obviously still single output discriminant
				d = setfeatlab(d,d.featlab(1:2,:));
			end             % needs conversion to two-classes before sigmoid
			d = sigm(d);
			d = normm(d);
		end
		%DXD finally, apply the cost matrix when it is defined in the
		%mapping/dataset:
		d = costm(d);
	end

elseif isa(a,'mapping') & isa(b,'dataset')
	error('Datasets should be given as first argument')

else
	error('Data type not supported')
end

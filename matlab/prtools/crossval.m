%CROSSVAL Error/performance estimation by cross validation (rotation)
% 
%   [ERR,CERR,NLAB_OUT] = CROSSVAL(A,CLASSF,N,1,TESTFUN)
%   [ERR,STDS]          = CROSSVAL(A,CLASSF,N,NREP,TESTFUN)
%   R                   = CROSSVAL(A,[],N,0)
%
% INPUT
%   A          Input dataset
%   CLASSF     The untrained classifier to be tested.
%   N          Number of dataset divisions (default: N==number of
%              samples, leave-one-out)
%   NREP       Number of repetitions (default: 1)
%   TESTFUN    Mapping,evaluation function (default classification error)
%
% OUTPUT
%   ERR        Average test error or performance weighted by class priors.
%   CERR       Unweighted test errors or performances per class
%   NLAB_OUT   Assigned numeric labels
%   STDS       Standard deviation over the repetitions.
%   R          Index array with rotation set
%
% DESCRIPTION
% Cross validation estimation of the error or performance (defined by TESTFUN)
% of the untrained classifier CLASSF using the dataset A. The set is randomly
% permutated and divided in N (almost) equally sized parts. The classifier
% is trained on N-1 parts and the remaining part is used for testing.  This
% is rotated over all parts. ERR is their weighted avarage over the class
% priors. CERR are the class error frequencies.  A and/or CLASSF may be
% cell arrays of datasets and classifiers. In that case ERR is an array
% with on position ERR(i,j) the error or performance of classifier j for
% dataset i. In this mode CERR and NLAB_OUT are returned in cell arrays.
%
% In case NREP > 1 the mean error(s) over the repetitions is returned in ERR
% and the standard deviations in the observed errors in STDS.
%
% In case NREP == 0 an index array is returned pointing to a fold for every
% object. No training or testing is done. This is useful for handling
% training and testing outside CROSSVAL.
% 
% See also DATASETS, MAPPINGS, TESTC

% Copyright: D.M.J. Tax, R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: crossval.m,v 1.10 2008/07/03 09:05:50 duin Exp $

function [err,cerr,nlabout] = crossval(data,classf,n,nrep,testf,fid)

	prtrace(mfilename);

  if nargin < 6, fid = []; end
  if nargin < 5, testf = []; end
	if nargin < 4, nrep = []; end
	if nargin < 3, n = []; end
  
  if ~ismapping(testf) & isempty(fid) % correct for old call without testf
    fid = testf; testf = [];
  end
	
	if iscell(data) % generate prior warnings now
		for j=1:length(data)
			data{j} = setprior(data{j},getprior(data{j}));
		end
	else
		data = setprior(data,getprior(data));
	end
	
	warnlevel = prwarning;
	prwarning(0);
	
	% datasets or classifiers are cell arrays
	if iscell(classf) | iscell(data)

		seed = rand('state');
		if ~iscell(classf), classf = {classf}; end
		if ~iscell(data), data = {data}; end
		if isdataset(classf{1}) & ismapping(data{1}) % correct for old order
			dd = data; data = classf; classf = dd;
		end
		numc = length(classf);
		numd = length(data);
		cerr = cell(numd,numc);
		nlab_out = cell(numd,numc);

		s1 = sprintf('crossval: %i classifiers: ',numc);
		prwaitbar(numc,s1);

		e = zeros(numd,numc);
    
		for jc = 1:numc

			prwaitbar(numc,jc,[s1 getname(classf{jc})]);
			s2 = sprintf('crossval: %i datasets: ',numd);
			prwaitbar(numd,s2);
			
 			for jd = 1:numd
				prwaitbar(numd,jd,[s2 getname(data{jd})]);
				rand('state',seed);
				[ee,cc,nn] = feval(mfilename,data{jd},classf{jc},n,nrep,testf);
				e(jd,jc) = ee;
				cerr(jd,jc) = {cc};
				nlabout(jd,jc) = {nn};
			end
			prwaitbar(0);

		end
		prwaitbar(0);
		if nrep > 1, cerr = cell2mat(cerr); nlabout = NaN; end
			
		if nargout == 0
			fprintf('\n  %i-fold cross validation result for',n);
			disperror(data,classf,e);
		end
		if nargout > 0, err  = e;  end

	else
		
		data = setprior(data,getprior(data)); % just to generate warning when needed
		
		if isempty(nrep), nrep = 1; end
		if nrep > 1
			
			s3 = sprintf('crossval: %i repetitions: ',nrep);
			prwaitbar(nrep,s3);
			ee = zeros(1,nrep);
			for j=1:nrep
				prwaitbar(nrep,j,[s3 int2str(j)]);
				[ee(j),ss,nlabout] = feval(mfilename,data,classf,n,1,testf);
			end
			prwaitbar(0);
			err = mean(ee);
			cerr = std(ee);
			nlabout = NaN;
			prwarning(warnlevel);
			return
		end

		if isdataset(classf) & ismapping(data) % correct for old order
			dd = data; data = classf; classf = dd;
		end
		isdataset(data);
		if nrep > 0, ismapping(classf); end
		[m,k,c] = getsize(data);
		lab = getlab(data);
		if isempty(n), n = m; end
		if n == m & ~isempty(testf)
			error('No external error routine allowed in case of leave-one-out cross validation')
		end

		if n > m
			warning('Number of batches too large: reset to leave-one-out')
			n = m;
		elseif n < 1
			error('Wrong size for number of cross-validation batches')
		end
		if (nrep > 0 & ~isuntrained(classf))
			error('Classifier should be untrained')
		end
		J = randperm(m);
		N = classsizes(data);

		% attempt to find an more equal distribution over the classes
		
		if all(N >= n)
			K = zeros(1,m);
			for i = 1:length(N)

				L = findnlab(data(J,:),i);

				M = mod(0:N(i)-1,n)+1;

				K(L) = M;

			end

		else
			K = mod(1:m,n)+1;

		end
		nlabout = zeros(m,1);
		
		if nrep == 0 % trick to return rotation set
			err = zeros(1,m);
			err(J) = K;
			prwarning(warnlevel);
			return
		end
		
    f = zeros(1,n);
		s4 = sprintf('crossval, %i-folds: ',n);
		prwaitbar(n,s4);
		
		for i = 1:n
			prwaitbar(n,i,[s4 int2str(i)]);
			OUT = find(K==i);
			JOUT=J(OUT);
			JIN = J; JIN(OUT) = [];
			train_data = data(JIN,:);
			%train_data = setprior(train_data,getprior(train_data));
			w = train_data*classf; % training
			                       % testing
			testres = data(JOUT,:)*w;
			if ~isempty(testf)
      	f(i) = testres*testf;
			end
			testout = testres*maxc;
			[mx,nlabout(JOUT)] = max(+testout,[],2);
                              % nlabout contains class assignments
		end

		prwaitbar(0);
		%correct for different order of lablist and labels assigned by
		%classifier. Assume this is constant over folds.
		nlist = renumlab(getfeatlab(testout),getlablist(data));
		nlabout = nlist(nlabout);
		if isempty(testf)
			e = zeros(1,c);
			for j=1:c
				J = findnlab(data,j);
				e(j) = sum(nlabout(J)~=j)/length(J);
    	end
			e = e*getprior(data,0)';
		else
			e = mean(f); % f already weighted by class priors inside testf
		end
		
		if nargout > 0
			err  = e;
			if isempty(testf)
				cerr = f;
			else
				cerr = [];
				nlabout = [];
			end
		else
			disp([num2str(n) '-fold cross validation error on ' num2str(size(data,1)) ' objects: ' num2str(e)])
		end

	end
	
	prwarning(warnlevel);
	return

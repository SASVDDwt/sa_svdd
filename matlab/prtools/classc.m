%CLASSC Convert mapping to classifier
%
%  W = CLASSC(W)
%  W = W*CLASSC
%
% INPUT
%  W  Any mapping or dataset
%	
% OUTPUT
%  W  Classifier mapping or normalized dataset: outputs/features sum to 1
%
% DESCRIPTION
% The mapping W is converted into a classifier by normalizing the outputs:
% the sum of the outputs for one sample equals 1. It is assumed that W 
% already generates densities or non-normalized posterior probability 
% estimates. This is true for neural networks or for two-class discriminants 
% calling CNORMC.
% A one-dimensional map is converted into a two-class classifier, provided
% that during the construction a class label was supplied. If not, the map
% cannot be converted and an error is generated.
% If D = A*W is the result of a mapping, it is normalized (the features of
% each sample sum to 1). This is identical to D = NORMM(D).
%
% SEE ALSO
% MAPPINGS, DATASETS

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: classc.m,v 1.4 2008/07/29 11:50:32 duin Exp $

function w = classc(w,flag)

	prtrace(mfilename);

	if nargin < 2, flag = 0; end % flag forces non=combiner behavior avoiding recursion
	
	if (nargin == 0)

		% Untrained mapping.
		w = mapping('classc','combiner',flag);			

	elseif (ismapping(w))

		% If mapping is stacked or parallel, recurse over the individual
		% sub-mappings and call CLASSC for each of them.
		if ((isstacked(w)) | (isparallel(w))) & (flag == 0)
			v = cell(1,length(w.data));
			for j = 1:length(w.data)
				if ismapping(w.data{j}) % the parallel combiner may have nonmapping data
					v{j} = feval(mfilename,w.data{j});
				else
					v{j} = w.data{j};
				end
			end
			w = setdata(w,v);
			w = feval(mfilename,w,1);  % and here CLASSC is called for the combiner avoiding recursion
		else
			conv = get(w,'out_conv');
			if (conv < 1)
				% Set the "normalization" bit in the mapping's output conversion flag
				w = set(w,'out_conv',conv+2);			
			else
				prwarning(3,'mapping is already a classifier');
			end;
		end
	elseif (isdataset(w))
		w = w*normm;
		w = w*costm;
	else
		error('input should be mapping or dataset');
	end

return

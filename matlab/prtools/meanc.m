%MEANC Mean combining classifier
% 
%   W = MEANC(V)
%   W = V*MEANC
%
% INPUT
%   V    Set of classifiers (optional)
%
% OUTPUT
%   W    Mean combiner
%
% DESCRIPTION
% If V = [V1,V2,V3, ... ] is a set of classifiers trained on the same
% classes and W is the mean combiner: it selects the class with the mean of
% the outputs of the input classifiers. This might also be used as
% A*[V1,V2,V3]*MEANC in which A is a dataset to be classified.
% 
% If it is desired to operate on posterior probabilities then the input
% classifiers should be extended like V = V*CLASSC;
%
% For affine mappings the coefficients may be averaged instead of the
% classifier results by using AVERAGEC.
% 
% The base classifiers may be combined in a stacked way (operating in the
% same feature space by V = [V1,V2,V3, ... ] or in a parallel way
% (operating in different feature spaces) by V = [V1;V2;V3; ... ]
%
% EXAMPLES
% PREX_COMBINING
%
% SEE ALSO
% MAPPINGS, DATASETS, VOTEC, MAXC, MINC, MEDIANC, PRODC,
% AVERAGEC, STACKED, PARALLEL

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands 

% $Id: meanc.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function w = meanc(p1)

	type = 'mean';               % define the operation processed by FIXEDCC.
	name = 'Mean combiner';      % define the name of the combiner.

	% this is the general procedure for all possible
	% calls of fixed combiners handled by FIXEDCC
	if nargin == 0
		w = mapping('fixedcc','combiner',{[],type,name});
	else
		w = fixedcc(p1,[],type,name);
	end

	if isa(w,'mapping')
		w = setname(w,name);
	end

	return

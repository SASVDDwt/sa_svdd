%PARZENML Optimum smoothing parameter in Parzen density estimation.
%         Soft label version
% 
% 	H = PARZENML(A,FID)
% 
% INPUT	
%   A    input dataset
%   FID  File ID to write progress to (default [], see PRPROGRESS)
%
% OUTPUT
%   H    scalar smoothing parameter
%
% DESCRIPTION
% Maximum likelihood estimation for the smoothing parameter H in the 
% Parzen denstity estimation of the data in A. A leave-one out 
% maximum likelihood estimation is used. 
%
% This routine does not use class information and computes a single
% smoothing parameter. It may be profitable to scale the data before
% calling it. eg. WS = SCALEM(A,'variance'); A = A*WS; If desired,
% remove unlabeled objects first, e.g. by SELDAT.
% 
% SEE ALSO
% DATASETS, MAPPINGS, SCALEM, SELDAT, PARZENML, PARZENDC, PRPROGRESS

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: parzenml.m,v 1.7 2004/01/21 20:45:44 bob Exp $

function h = parzenml(A,fid)

	prtrace(mfilename);
	
	if nargin < 2, fid = []; end
	
	SS = gettargets(setlabtype(A,'soft'));
	[m,k,c] = getsize(A);
	DD= distm(+A) + diag(1e70*ones(1,m));
	E = min(DD);
	h = zeros(c,1);
	h0 = sqrt(max(E));    % initial estimate of h
	
	for j=1:c
		
		S = SS(:,j);
		h1 = h0;
		F1 = derl(DD,E,h1,k,S); % derivative

		prprogress(fid,'\nparzenml: ML Smoothing Parameter Optimization\n')
		prprogress(fid,'  h = %5.3f   F = %8.3e \n',h1,F1);
		if abs(F1) < 1e-70 
			h = h1;
			prwarning(4,'jump out\n');
			return;
		end
	
		a1 = (F1+m*k)*h1*h1;
		h2 = sqrt(a1/(m*k));  % second guess
		F2 = derl(DD,E,h2,k,S); % derivative

		prprogress(fid,'  h = %5.3f   F = %8.3e \n',h2,F2);
		if (abs(F2) < 1e-70) | (abs(1e0-h1/h2) < 1e-6) 
			h(j) = h2;
			prwarning(4,'jump out\n');
			break;
		end
	
		% find zero-point of derivative to optimize h^2
		% stop if improvement is small, or h does not change significantly
	
		alf = 1;
		while abs(1e0-F2/F1) > 1e-4 & abs(1e0-h2/h1) > 1e-3 & abs(F2) > 1e-70
			h3 = (h1*h1*h2*h2)*(F2-F1)/(F2*h2*h2-F1*h1*h1);
			if h3 < 0 % this should not happen
				h3 = sqrt((F2+m*k)*h2*h2/(m*k));
			else
				h3 = sqrt(h3);
			end
			h3 = h2 +alf*(h3-h2);
			F3 = derl(DD,E,h3,k,S);
			prprogress(fid,'  h = %5.3f   F = %8.3e \n',h3,F3);
			F1 = F2; F2 = F3;
			h1 = h2; h2 = h3;
			alf = alf*0.99; % decrease step size
		end
		h(j) = h2;
		prprogress(fid,'parzenml finished')
	end

return

function F = derl(DD,E,h,k,S)
	% computation of the likelihood derivative for Parzen density
	% given distances D and their object minima E (for increased accuracy)
	% S are the object weigths
	c = size(S,2);                  % number of classes
	m = size(DD,1);
	Y = (DD-repmat(E,m,1))/(2*h*h); % correct for minimum distance to save accuracy
	IY = find(Y<20);                % take small distance only, others don't contribute
	F = 0;
	for j=1:c
		P = zeros(m,m);
		P(IY) = exp(-Y(IY));
		PP = S(:,j)'*P';
		FU = repmat(realmax,1,m);
		J = find(PP~=0);  
		FU(J) = S(J,j)'./PP(J);
		K = find(S(:,j)==0);
		FU(K) = zeros(1,length(K));
		FF = (DD.*P)*S(:,j);
		F = F + (FU*FF)./(h*h);
	end
	F = F - sum(S(:))*k;
return

%TESTD Replaced by TESTC

function [errors,class_errors] = testd(a,w)

global TESTD_REPLACED_BY_TESTC

if isempty(TESTD_REPLACED_BY_TESTC)
	disp([newline 'TESTD has been replaced by TESTC, please use it'])
	TESTD_REPLACED_BY_TESTC = 1;
end

if nargin == 0
	errors = testc;
elseif nargin == 1
	[errors,class_errors] = testc(a);
else
	[errors,class_errors] = testc(a,w);
end

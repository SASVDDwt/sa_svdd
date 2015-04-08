%IM_FILL_NORM Fill and normalize image for display puproses
%
%   B = IM_FILL_NORM(A,N)
%
%Low level routine for the DATAFILE/SHOW command to display non-square
%images of the datafile A, inside square of NxN pixels. Empty areas are
%filled with gray.

function b = im_fill_norm(a,n)

if isa(a,'dataset')
	isobjim(a);
	outsize = [n n getfeatsize(a,3)];
  b = filtim(a,'im_fill_norm',{n},outsize);
else
	a = double(a);
	[x,y,p] = size(a);
	mx = max(a(:));
	mn = min(a(:));
	b = 0.5*ones(n,n,p);
	%b = ones(n,n,p);
	z = floor(n/2);
	if x > y
		a = imresize(a,round(n*[x,y]/x));
		k = size(a,2);
		s = floor((n-k)/2)+1;
		b(:,s:s+k-1,:) = a;
	else
		a = imresize(a,round(n*[x,y]/y));
		k = size(a,1);
		s = floor((n-k)/2)+1;
		b(s:s+k-1,:,:) = a;
	end
	a = (a - mn)/(mx-mn+eps);
end

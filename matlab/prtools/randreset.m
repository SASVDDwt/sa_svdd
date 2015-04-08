
function output = randreset(state)

if nargin < 1
   state = 1;
end

randstate = cell(1,2);
randstate{1} = rand('state');
randstate{2} = randn('state');

if iscell(state)
   rand('state',state{1});
   randn('state',state{2});
else
   rand('state',state);
   randn('state',state);
end

if nargout > 0
  output = randstate;
end
   

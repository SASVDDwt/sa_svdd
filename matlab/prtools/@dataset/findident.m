%FINDIDENT Determine indices of objects having specified identifiers
%
%	   J = FINDIDENT(A,IDENT,FIELD)
%
% INPUT
%   A      Dataset
%   IDENT  Object identifiers, see SETIDENT
%   FIELD  Desired field, default 'IDENT'.
%
% If IDENT is a single object identifier then J is a vector with indices
% to all objects with that identifier in the specified FIELD.
%
% If IDENT is a set of object identifiers, then J is a column vector of
% cells, such that J{n} contains a vector with indices to all objects having
% identifier IDENT{n}.

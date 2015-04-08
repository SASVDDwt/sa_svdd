%DATASET Dataset class constructor
%
%    A = DATASET(DATA,LABELS,VARARGIN)
%
% A DATASET object is constructed from:
%
% DATA    size [M,K], a set of M datavectors of K features
%                     a cell array of datasets will be concatenated. 
% LABELS  size [M,N]  array with labels for the M datavectors.
%                     LABELS should be either integers or character strings.
%                     Choose single characters for the fastest implementation.
%                     Numeric labels with value NaN or character labels
%                     with value CHAR(0) are interpreted as missing labels.
%                     See also RENUMLAB.
%
% Other parameter fields may be set by
%
%    A = DATASET(DATA,LABELS,'field1',VALUE1,'field2',VALUE2, ...)
%
% The following parameter fields are possible:
%
% FEATLAB size [K,F]  array with labels for the K features
% FEATDOM size [K]    cell array with domain description for the K features
% PRIOR   size [C,1]  prior probabilities for each of the C classes
%                     PRIOR = 0: all classes have equal probability 1/C
%                     PRIOR = []: all datavectors are equally probable
% COST    size [C,C+1] Classification cost matrix. COST(I,J) are the costs
%                     of classifying an object from class I as class J.
%                     Column C+1 generates an alternative reject class and
%                     may be omitted, yielding a size of [C,C]. 
%                     An empty cost matrix, COST = [] (default) is interpreted
%                     as COST = ONES(C) - EYE(C) (identical costs of
%                     misclassification).
% LABLIST size [C,N]  class labels corresponding to the unique labels found
%                     in LABELS and thereby to the classes in the dataset.
%                     The order of the items in LABLIST corresponds to the
%                     apriori probablities stored in PRIOR. LABLIST should
%                     only be given explicitely if PRIOR is given and if it
%                     is not equal to 0 and not empty.
% LABTYPE             String defining the label type,
%                     'crisp' for defining classes by integers or strings
%                     'soft' for defining memberships to classes. In this
%                             case LABELS should be a MxC array with numbers
%                             between 0 and 1.
%                     'targets' for defining regression type target values.
%                             Labels should be a MxN numeric array for
%                             defining N targets per object.
% OBJSIZE             number of objects, or vector with its shape. This is
%                     useful if the set of objects can be interpreted as an
%                     image (objects are pixels).
% FEATSIZE            number of features, or vector with its shape. This is
%                     useful if the set of features can be interpreted as an
%                     image (features are pixels).
% IDENT  [M,1]        Cell array, identifier for objects. 
% NAME                String with dataset name
% USER                User definable variable
%
% These parameters are parsed and stored in the following fields:
%
% A.DATA    = data
% A.NLAB    = numeric labels, index in lablist
% A.FEATLAB = feature labels
% A.FEATDOM = feature domains
% A.PRIOR   = prior probabilities
% A.COST    = classification cost matrix
% A.LABLIST = labels of the classes
% A.TARGETS = dataset with soft labels or targets
% A.LABTYPE = label type: 'crisp','soft' or 'target'
% A.OBJSIZE = number of objects or vector with its shape
% A.FEATSIZE= number of features or vector with its shape
% A.IDENT   = identifier for objects (integer)
% A.VERSION = PRTools version used for creating dataset
% A.NAME    = string with name of the dataset
% A.USER    = user field
%
% Objects of type MEASUREMENT or old DATASET definitions given by
% by a structure can be converted by the DATASET constructor.
%
% Data can be added or changed in an existing dataset by:
% SET, SETDATA, SETFEATLAB, SETFEATDOM, SETFEATSIZE, SETIDENT,
% SETLABELS, SETLABLIST, SETLABTYPE, SETNAME, SETNLAB, SETOBJSIZE,
% SETPRIOR, SETCOST, SETTARGETS, SETUSER.
%
% Data can be retrieved from a dataset by:
% GET, GETDATA, GETFEATLAB, GETFEATDOM, GETFEATSIZE, GETIDENT, GETLABELS,
% GETLABLIST, GETLABTYPE, GETNAME, GETNLAB, GETOBJSIZE, GETPRIOR, GETCOST,
% GETSIZE, GETTARGETS, GETUSER, GETVERSION, FINDIDENT, FINDLABELS,
% FINDNLAB.
%
% Shortcuts for retrieving the datafield A.DATA are DOUBLE(A) and +A.
%
% SEE ALSO DATASETS

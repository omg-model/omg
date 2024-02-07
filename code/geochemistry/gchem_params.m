%% GEOCHEMISTRY FUNCTIONS
gchem_fcns                                =    geochemistry_functions; 


%% DEFAULT GEOCHEMISTRY PARAMETERS

%% FOR CO2SYS

% Set pH scale
%  1 = Total scale
%  2 = Seawater scale
%  3 = Free scale  
%  4 = NBS scale
gchem_pars.pHscale          = 2;

%   K1K2CONSTANTS     : scalar or vector of size n (***)
%  (***) Each element must be an integer,
%        indicating the K1 K2 dissociation constants that are to be used:
%   1 = Roy, 1993											T:    0-45  S:  5-45. Total scale. Artificial seawater.
%   2 = Goyet & Poisson										T:   -1-40  S: 10-50. Seaw. scale. Artificial seawater.
%   3 = HANSSON              refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   4 = MEHRBACH             refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   5 = HANSSON and MEHRBACH refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   6 = GEOSECS (i.e., original Mehrbach)					T:    2-35  S: 19-43. NBS scale.   Real seawater.
%   7 = Peng	(i.e., originam Mehrbach but without XXX)	T:    2-35  S: 19-43. NBS scale.   Real seawater.
%   8 = Millero, 1979, FOR PURE WATER ONLY (i.e., Sal=0)	T:    0-50  S:     0.
%   9 = Cai and Wang, 1998									T:    2-35  S:  0-49. NBS scale.   Real and artificial seawater.
%  10 = Lueker et al, 2000									T:    2-35  S: 19-43. Total scale. Real seawater.
%  11 = Mojica Prieto and Millero, 2002.					T:    0-45  S:  5-42. Seaw. scale. Real seawater
%  12 = Millero et al, 2002									T: -1.6-35  S: 34-37. Seaw. scale. Field measurements.
%  13 = Millero et al, 2006									T:    0-50  S:  1-50. Seaw. scale. Real seawater.
%  14 = Millero et al, 2010									T:    0-50  S:  1-50. Seaw. scale. Real seawater.
gchem_pars.K1K2             = 4;

%  (****) Each element must be an integer that
%         indicates the KSO4 dissociation constants that are to be used,
%         in combination with the formulation of the borate-to-salinity ratio to be used.
%         Having both these choices in a single argument is somewhat awkward,
%         but it maintains syntax compatibility with the previous version.
%  1 = KSO4 of Dickson & TB of Uppstrom 1979  (PREFERRED)
%  2 = KSO4 of Khoo    & TB of Uppstrom 1979
%  3 = KSO4 of Dickson & TB of Lee 2010
%  4 = KSO4 of Khoo    & TB of Lee 2010
gchem_pars.KSO4             = 1;

%%
gchem_pars.constants_method      = 'OMG';
gchem_pars.interpolate_constants = true;

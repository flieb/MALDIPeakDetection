% ALLDIV.M
%  PROPDIV   - Proper divisors of integer
% 
%  Usage:      d=propdiv(n)
% 
%  Input:      n  integer
%  Output:     d  vector of proper divisors of n
% 
% ---------------------------------------------------------------
%
% COPYRIGHT : (c) NUHAG, Dept.Math., University of Vienna, AUSTRIA
%             http://nuhag.eu/
%             Permission is granted to modify and re-distribute this
%             code in any manner as long as this notice is preserved.
%             All standard disclaimers apply.
%
function d=propdiv(n)
% PROPDIV   - Proper divisors of integer
%
% Usage:      d=propdiv(n)
%
% Input:      n  integer
% Output:     d  vector of proper divisors of n

% N. Kaiblinger, 2000

d=n./(n-1:-1:2);
d=d(d==round(d));

% end of file

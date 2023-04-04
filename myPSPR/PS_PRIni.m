%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%  Copyright   Xin-Guang Zhu, Yu Wang, Donald R. ORT and Stephen P. LONG
%  CAS-MPG Partner Institute for Computational Biology, SIBS/CAS, China 
%  Institute of Genomic Biology and Department of Plant Biology,
%  University of Illinois at Urbana Champaign, United States
%  This file is part of e-photosynthesis. All Rights Reserved.
 
%   e-photosynthesis is distributed for academic research only. 
%   For commercial purpose, please apply for different liscence. 
%   By using this software, you are automatically accepting the academic
%   free licence http://opensource.org/licenses/afl-3.0.php
%   See the content of the licence in the word document: 
%   Academic Free Licence.doc
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function PS_PRs = PS_PRIni(params)

PSs = PSInitial(params);
PrS = PRinitial(params);

PS_PRs = zeros(24,1);

for m=1:4
    PS_PRs(m) = PSs(m);
end

for m = 5:14
    PS_PRs(m) = PSs(m+1);
end

for m = 15:16
    PS_PRs(m) = PrS(m-14);
end

for m = 17:23
    PS_PRs(m) = PrS(m-13);
end

PS_PRs(24) = PSs(5);
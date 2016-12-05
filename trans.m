## Copyright (C) 2016 Gonzalo Rodríguez Prieto <gonzalo.rprieto@uclm.es>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.
###########################################################################################
## -*- texinfo -*-
## @deftypefn {Function File} {} trans ()
##  Transforms text chains in sclaing factors, following the S.I. system
##
## @code{trans (@var{multiplier})} 
##  returns the nucmerical scaling for the string @var{multiplier}.
##  If not @var{output} is given, the stdout is used.
##  No empty strings are allowed.
## @end deftypefn

## Made by: Gonzalo Rodríguez Prieto (gonzalo.rprietoATuclm.es)

function output = trans(multiplier)

if length(multiplier) == 0
 error("function trans: No empty strings are allowed.");
endif; 


switch (multiplier)
 case "mm"
	output = 1e-3;
 case "S.I."
    output = 1; %No changes needed, already SI. units and correctly scaled...
 case "um"
    output = 1e-6;
 case "k"
    output = 1e3;
 case "g/cm3"
    output = 1e3;
 case "g/mol"
    output = 1e-3;
 case "K"
    output = 1;
 case "J"
    output = 1;           
 otherwise
   output = 1;
endswitch;
 
endfunction;

%And that's that's all folks!!!

function [T_load] = findTload(h5)
%%% This program uses the ""Function"" command to call inputs from theta 2 and h5 and 
%%% runs it into the program to obatin the value of h5

if h5>0
    T_load=-5; %ft*lbf constant load
elseif h5<0
    T_load=5;  %ft*lbf constant load
end
end
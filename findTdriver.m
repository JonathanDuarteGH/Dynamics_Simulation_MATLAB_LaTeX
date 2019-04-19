function [T_driver] = findTdriver(theta2dot, T_stall, w_max)
%%% This program uses the ""Function"" command to call inputs from 
%%% theta2dot, T_stall and 
%%% w_max and runs it into the program to obatin the value of h5
T_driver=T_stall.*(1-theta2dot/w_max);
end
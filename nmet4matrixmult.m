function [t3r,t4r,t5r]= nmet4matrixmult(r1,r2,r3,r4,t2r,t3r,t4r,t5r,t5ri,t4radi,t3radi)
%%% This program uses the ""Function"" command to call inputs from x_o and 
%%% runs it into the program until the error is less than or equal to 
%%% the tolerance.
%%% This program was used in Programming Problem 2; Its now being implemented in
%%% Design Set 2. 
%%% r1, r2, r3, r4 are all known link lengths (in). 
%%% t2 will be in a range of 0-2pi in terms of radians.

x_o = [t3r;t4r;t5r]; %;t5r] % defining the vector of inital guess which 
                 % will already in radians once the function 
                 % is called in the main program
% setting tolerance for convergence
tol= .1;
% setting initial error greater than tolerance to begin loop
error= 10;
while error>tol
    % setting sines and cosines of angles t2, t3, t4
    st2= sin(t2r);
    ct2= cos(t2r);
    st3= sin(t3r);
    ct3= cos(t3r);
    st4= sin(t4r);
    ct4= cos(t4r);
    f1= r2*ct2+r3*ct3-r4*ct4-r1;   
    f2= r2*st2+r3*st3-r4*st4;
    f3= (t5r-t5ri)-2*(t4r-t4radi)+(t3r-t3radi);
    f= [f1; f2; f3];
    dfdt3= [-r3*st3; r3*ct3; 1];  
    dfdt4= [r4*st4; -r4*ct4; -2]; 
    dfdt5= [0; 0; 1];
    A= [dfdt3,dfdt4,dfdt5];
    x= -pinv(A)*f+ x_o;
    t3r= x(1,1);
    t4r= x(2,1);
    t5r= x(3,1);
    %t5r= 2*(t4r - t4radi) - (t3r - t3radi) + t5ri;
    error = norm(x-x_o);
    x_o=x;
end
end

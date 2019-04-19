function [f_g5xp, f_g5yp, h5p, h4p, h3p, f_g5y, f_g5x, h5, h4, h3, r_g5y, r_g5x] =...
    first_second_kc(r5, theta5, r2, r3, r4, t4r, t3r, t2r, count)

% Half of WIPER BLADE Length
    r_g5=r5/2;
    r_g5x=r_g5*cos(theta5);
    r_g5y=r_g5*sin(theta5);
    %%%% 1st Order KC
    h3 = r2.*sin(t4r-t2r)./(r3.*sin(t3r-t4r));
    h4 = r2.*sin(t3r-t2r)./(r4.*sin(t3r-t4r));
    h5 = 2*h4(count)-h3(count);
    %%%% 1st Order KC WIPER BLADE
    f_g5x=-1*r_g5*sin(theta5).*h5;
    f_g5y=r_g5*cos(theta5).*h5;
    %%%% 2nd Order KC
    J11=-1*r3*sin(t3r); J12=r4*sin(t4r); J13=0;
    J21=r3*cos(t3r); J22=-1*r4*cos(t4r); J23=0;
    J31=1; J32=2; J33=1;
    J = [J11 J12 J13; J21 J22 J23; J31 J32 J33];
    A11=r2*cos(t2r)+(r3*cos(t3r)*h3.^2)-(r4*cos(t4r)*h4.^2);
    A12=r4*sin(t4r);
    A13=0;
    A21=r2*sin(t2r)+(r3*sin(t3r)*h3.^2)-(r4*sin(t4r)*h4.^2);
    A22=-1*r4*cos(t4r);
    A23=0;
    A31=0; A32=-2; A33=1;
    A = [A11 A12 A13; A21 A22 A23; A31 A32 A33];
    B11=-1*r3*sin(t3r);
    B12=r2*cos(t2r)+(r3*cos(t3r)*h3.^2)-(r4*cos(t4r)*h4.^2);
    B13=0;
    B21=r3*cos(t3r);
    B22=r2*sin(t2r)+(r3*sin(t3r)*h3.^2)-(r4*sin(t4r)*h4.^2);
    B23=0;
    B31=1; B32=0; B33=1;
    B = [B11 B12 B13; B21 B22 B23; B31 B32 B33];
    C11=-1*r3*sin(t3r);
    C12=r4*sin(t4r);
    C13=r2*cos(t2r)+(r3*cos(t3r)*h3.^2)-(r4*cos(t4r)*h4.^2);
    C21=r3*cos(t3r);
    C22=-1*r4*cos(t4r);
    C23=r2*sin(t2r)+(r3*sin(t3r)*h3.^2)-(r4*sin(t4r)*h4.^2);
    C31=1; C32=-2; C33=0;
    C = [C11 C12 C13; C21 C22 C23; C31 C32 C33];
    h3p = det(A)/det(J);
    h4p = det(B)/det(J);
    h5p = det(C)/det(J);
    %%%% 2nd Order KC WIPER BLADE
    f_g5xp=-1*r_g5*cos(theta5).*(h5.^2)-r_g5*sin(theta5).*h5p;
    f_g5yp=-1*r_g5*sin(theta5).*(h5.^2)+r_g5*cos(theta5).*h5p;
    
end
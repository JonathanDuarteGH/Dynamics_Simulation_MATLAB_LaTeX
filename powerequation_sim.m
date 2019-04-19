function [Sum_CAP_A5, Sum_CAP_B5] =...
    powerequation_sim(h2, m_wb_slug, f_g5x, theta2, f_g5y, I_gwb, f_g5xp,...
    f_g5yp, h5, h5p, I_gmw, I_gwg, I_gww)

    %%%% Defining the Orientation of the Motor Windings Worm Gear and Worm
    %%%% Wheel
    h_wheel=h2; h_wormg=83*h_wheel; h_motor=h_wormg;
    %%%% A and B for WIPER BLADE %%Slugs ft^{2}%%
    CAP_A1=m_wb_slug.*(f_g5x.^2+f_g5y.^2)+I_gwb.*h5.^2;
    CAP_B1=m_wb_slug.*(f_g5x.*f_g5xp+f_g5y.*f_g5yp)+I_gwb.*h5.*h5p;
    %%%% A and B for Motor Windings %%slug ft^{2}%%
    CAP_A2(length(theta2))=I_gmw.*h_motor.^2;
    CAP_B2(length(theta2))=I_gmw.*h_motor.*0;
    %%%% A and B for Worm Gear %%slug ft^{2}%%
    CAP_A3(length(theta2))=I_gwg.*h_wormg.^2;
    CAP_B3(length(theta2))=I_gwg.*h_wormg.*0;
    %%%% A and B for Worm Wheel %%slug ft^{2}%%
    CAP_A4(length(theta2))=I_gww.*h_wheel.^2;
    CAP_B4(length(theta2))=I_gww.*h_wheel.*0;
    %%%% Sum of the A's and B's %%slug ft^{2}%%
    Sum_CAP_A5=CAP_A1+CAP_A2+CAP_A3+CAP_A4;
    Sum_CAP_B5=CAP_B1+CAP_B2+CAP_B3+CAP_B4;
    
end
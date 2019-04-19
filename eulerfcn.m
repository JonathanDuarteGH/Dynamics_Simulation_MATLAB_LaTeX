function [theta5dot, deltatheta2, dt, deltatheta2dot, dsi] =...
    eulerfcn(theta2_intial, theta2dot, theta2dotdot, h5)
dsi = theta2_intial*(pi/180);

dt = (-theta2dot+sqrt(theta2dot^2+(2*theta2dotdot*dsi)))/theta2dotdot;

deltatheta2dot = theta2dotdot*dt;
deltatheta2 = (theta2dot*dt)+(0.5*theta2dotdot*(dt^2));

theta5dot = h5*theta2dot;
end
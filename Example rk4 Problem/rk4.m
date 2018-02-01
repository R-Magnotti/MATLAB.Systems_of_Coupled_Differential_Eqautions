%Subroutine for rkTestAnharm
%This is the Runge-Kutta Stepper
function[T,Y,Z] = rk4(f,g,a,b,ya,za,m)



h = (b-a)/m; %step size



T = zeros(1,m+1); %time vector (stores all t values)

Y = zeros(1,m+1); %position vector (stores all position y values)

Z = zeros(1,m+1); %velcotiy vector (stores all velocity v values)



T(1) = a; %initial time

Y(1) = ya; %intial position

Z(1) = za; %initial velocity



for j=1:m, %iterates one time step (entails 4 approximations)

    
    tj = T(j); %increase j value per iteration of loop

    yj = Y(j);

    zj = Z(j);

    

    k1 = h*feval(f,tj,yj,zj); %h(run)*(slope of position (y) tangent line) = approx. rise from last iteration

    l1 = h*feval(g,tj,yj,zj); %h(run)*(slope of velocity (z) tangent line) = approx. rise from last iteration

    k2 = h*feval(f,tj+(h/2),yj+(k1/2), zj+(l1/2)); % "" but using updated values of y and z (half steps)

    l2 = h*feval(g,tj+(h/2),yj+(k1/2),zj+(l1/2)); % "" (half steps)

    k3 = h*feval(f,tj+(h/2),yj+(k2/2),zj+(l2/2)); % "" (half steps)

    l3 = h*feval(g,tj+(h/2),yj+(k2/2),zj+(l2/2)); % "" (half steps)

    k4 = h*feval(f,tj+h,yj+k3,zj+l3); % ""

    l4 = h*feval(g,tj+h,yj+k3,zj+l3); % ""

    

    Y(j+1) = yj + (k1 + 2*k2 + 2*k3 + k4)/6; %computes new y value by averaging all 4 approximations of y

    Z(j+1) = zj + (l1 + 2*l2 + 2*l3 + l4)/6; %"" z

    T(j+1) = a + h*j; %updates by taking another step

end


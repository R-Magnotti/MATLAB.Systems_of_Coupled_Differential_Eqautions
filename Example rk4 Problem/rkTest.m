%Anharmonic practice problem
%following eq.'s are how second order 
%ODE breaks down into two, coupled first order
%ODE's

k = 1.0; %spring constant

c = 0; %some constant

f = @(t,y,z) z; 

g = @(t,y,z) -k*y - c*y^3;



a=0; %intiial time

b=1000; %terminal time

ya=3; %initial position

za=0; %initial velocity

m=1500; %number of steps



[T,Y,Z] = rk4(f,g,a,b,ya,za,m); %returns T,Y,Z vectors 
                                %with each respective t,y,z values 
                                %at each time step


%plot(T,Y)

%Subroutine for AllSpinLattice Code
%This is the Runge-Kutta Stepper
function[T,S,E,HEFF] = TSLrk4(fMat,a,b,spinMat,m,J,heffMat,s)

% Basic Variables

h = (b-a)/m; %step size

T = zeros(1,m+1); %time vector (stores all t values)...(1xm+1) matrix

T(1) = a; %initial time

% 3D Spin matrix (spin xyz values -> rows, spin# -> columns, spin at t -> page)
% pages change w/ time steps

S = zeros(3,m+1,s);

for j = 1:s
    
    S(:,j,1) = spinMat(:,j); %storing intial spin values
    
end

% 3D heff matrix (heff xyz values -> rows, heff# -> columns, heff @ t -> page)
% pages change w/ time steps

HEFF = zeros(3,m+1,s);

for j = 1:s
    
    HEFF(:,j,1) = heffMat(:,j); %storing intial heff value per particle
    
end

% Vector with one row and m+1 columns to store the total energy of the 
% system per time step

E = zeros(1,m+1); %to store all Hamiltonian values

% matrix to represent the dot products of each matrix with two nearest
% neighbors (spin dot prod. value -> row, spin nearest neighbor pair # ->
% column)

sDotMat = zeros(1,m+1);
sDotMat(1,1) = dot(s1a,s2a); %dot prod. of initial s values    

E(:,1) = -J*sDot; %Initial energy
    
---------------------------------------------------------------------------
for j=1:m, %iterates one time step (entails 4 approximations)

           %like for loop in Java (for(j = 1; j = m; j++))

    

    tj = T(j); %increase j value per iteration of loop

    s1j = S1(:,j);
    s2j = S2(:,j);
    s3j = S3(:,j);
    
    f1 = fMat(1,j);
    f2
    
    heff1 = HEFFs1(:,j);
    heff2 = HEFFs2(:,j);    
    heff3 = HEFFs3(:,j);
    
    
    %start rk4 portion
    k1 = h*feval(f1,tj,heff1,s1j); %h(run)*(slope of spin1 (s1) tangent line) = approx. rise from last iteration
    l1 = h*feval(f2,tj,heff2,s2j); %h(run)*(slope of spin2 (s2) tangent line) = approx. rise from last iteration    
    m1 = h*feval(f3,tj,heff3,s3j);

    k2 = h*feval(f1,tj+(h/2),heff1+(l1/2),s1j+(k1/2)); % "" but using updated values of s1 and 2 (half steps)
    l2 = h*feval(f2,tj+(h/2),heff2+(k1/2),s2j+(l1/2)); % "" (half steps)
    m2 = h*feval(f3,tj+(h/2),heff3+(m1/2),s3j+(m1/2));

    k3 = h*feval(f1,tj+(h/2),heff1+(l2/2),s1j+(k2/2)); % "" (half steps)
    l3 = h*feval(f2,tj+(h/2),heff2+(k2/2),s2j+(l2/2)); % "" (half steps)
    m3 = h*feval(f3,tj+(h/2),heff3+(m2/2),s3j+(m2/2));

    k4 = h*feval(f1,tj+h,heff1+l3,s1j+k3); % ""
    l4 = h*feval(f2,tj+h,heff2+k3,s2j+l3); % ""
    m4 = h*feval(f3,tj+h,heff3+m3,s3j+m3);
    
    
    %takes old value of S1 and adds new value to it to update it
    s1j = s1j + (k1 + 2*k2 + 2*k3 + k4)/6; %computes new y value by averaging all 4 approximations of y
    s1j = normc(s1j);
    S1(:,j+1) = s1j;
    
    s2j = s2j + (l1 + 2*l2 + 2*l3 + l4)/6; %"" z
    s2j = normc(s2j);
    S2(:,j+1) = s2j;
    
    s3j = s3j + (m1 + 2*m2 + 2*m3 + m4)/6; %"" z
    s3j = normc(s3j);
    S3(:,j+1) = s3j;
    
    T(j+1) = a + h*j; %updates by taking another step

    HEFFs1(:,j+1) = -J*(s2j+s3j);
    HEFFs2(:,j+1) = -J*(s1j+s3j);      
    HEFFs3(:,j+1) = -J*(s1j+s2j);
    
    sDot1 = dot(s2j,s3j);
    sDot2 = dot(s1j,s3j);
    sDot3 = dot(s1j,s2j);
    
    sDotTot
    
    E(:,j+1) = -J*(sDot1 + sDot2 + sDot3); %Hamiltonian
end


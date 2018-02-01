%Subroutine for ThreeSpinLattice Code
%This is the Runge-Kutta Stepper
function[T,S1,S2,S3,E,HEFFs1,HEFFs2,HEFFs3] = ThSLrk4(f1,f2,f3,a,b,s1a,s2a,s3a,m,J,heff1,heff2,heff3)

h = (b-a)/m; %step size


T = zeros(1,m+1); %time vector (stores all t values)...(1xm+1) matrix

S1 = zeros(3,m+1); %Spin1 vector (stores all Spin1 values)
S2 = zeros(3,m+1); %Spin2 vector 
S3 = zeros(3,m+1);

S1(:,1) = s1a; %intial spin1 (particle1) value
S2(:,1) = s2a; 
S3(:,1) = s3a;

HEFFs1 = zeros(3,m+1); %h effective vector
HEFFs2 = zeros(3,m+1);
HEFFs3 = zeros(3,m+1);

HEFFs1(:,1) = heff1; %initial effective h of particle 1
HEFFs2(:,1) = heff2;
HEFFs3(:,1) = heff3;

T(1) = a; %initial time

E = zeros(1,m+1); %to store all Hamiltonian values

sDot1 = dot(s2a,s3a);
sDot2 = dot(s1a,s3a);
sDot3 = dot(s1a,s2a);
    
E(:,1) = -J*(sDot1+sDot2+sDot3); %Initial energy
    
for j=1:m, %iterates one time step (entails 4 approximations)

           %like for loop in Java (for(j = 1; j = m; j++))

    

    tj = T(j); %increase j value per iteration of loop

    s1j = S1(:,j);

    s2j = S2(:,j);
    
    s3j = S3(:,j);
    
    
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
    
    T(j+1) = a + h*j; %Updates by taking another step

    HEFFs1(:,j+1) = -J*(s2j+s3j);
    HEFFs2(:,j+1) = -J*(s1j+s3j);      
    HEFFs3(:,j+1) = -J*(s1j+s2j);
    
    sDot1 = dot(s2j,s3j);
    sDot2 = dot(s1j,s3j);
    sDot3 = dot(s1j,s2j);
    
    E(:,j+1) = -J*(sDot1 + sDot2 + sDot3); %Energy Hamiltonian
end


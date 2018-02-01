%Subroutine for TwoSpinLattice Code
%This is the Runge-Kutta Stepper
function[T,S1,S2,E,HEFFs1,HEFFs2] = TSLrk4(f,g,a,b,s1a,s2a,m,J, heffs1, heffs2)


h = (b-a)/m; %step size


T = zeros(1,m+1); %time vector (stores all t values)...(1xm+1) matrix

S1 = zeros(3,m+1); %Spin1 vector (stores all Spin1 values)
S2 = zeros(3,m+1); %Spin2 vector 

S1(:,1) = s1a; %intial spin1 (particle1) value
S2(:,1) = s2a; 

HEFFs1 = zeros(3,m+1); %h effective vector
HEFFs2 = zeros(3,m+1);

HEFFs1(:,1) = heffs1; %initial effective h of particle 1
HEFFs2(:,1) = heffs2;


T(1) = a; %initial time

E = zeros(1,m+1); %to store all Hamiltonian values

sDot = dot(s1a,s2a);
    
E(:,1) = -J*sDot; %Initial energy
    
for j=1:m, %iterates one time step (entails 4 approximations)

           %like for loop in Java (for(j = 1; j = m; j++))

    

    tj = T(j); %increase j value per iteration of loop

    s1j = S1(:,j);

    s2j = S2(:,j);
    
    heff1 = HEFFs1(:,j);

    heff2 = HEFFs2(:,j);
    
    k1 = h*feval(f,tj,heff1,s1j);
    %Original: k1 = h*feval(f,tj,s1j,s2j); %h(run)*(slope of spin1 (s1) tangent line) = approx. rise from last iteration

    l1 = h*feval(g,tj,heff2,s2j);
    %Original: l1 = h*feval(g,tj,s1j,s2j); %h(run)*(slope of spin2 (s2) tangent line) = approx. rise from last iteration

    k2 = h*feval(f,tj+(h/2),heff1+(l1/2),s1j+(k1/2)); % "" but using updated values of s1 and 2 (half steps)

    l2 = h*feval(g,tj+(h/2),heff2+(k1/2),s2j+(l1/2)); % "" (half steps)

    k3 = h*feval(f,tj+(h/2),heff1+(l2/2),s1j+(k2/2)); % "" (half steps)

    l3 = h*feval(g,tj+(h/2),heff2+(k2/2),s2j+(l2/2)); % "" (half steps)

    k4 = h*feval(f,tj+h,heff1+l3,s1j+k3); % ""

    l4 = h*feval(g,tj+h,heff2+k3,s2j+l3); % ""

    
    %takes old value of S1 and adds new value to it to update it
    s1j = s1j + (k1 + 2*k2 + 2*k3 + k4)/6; %computes new y value by averaging all 4 approximations of y
    s1j = normc(s1j);
    S1(:,j+1) = s1j;
    
    s2j = s2j + (l1 + 2*l2 + 2*l3 + l4)/6; %"" z
    s2j = normc(s2j);
    S2(:,j+1) = s2j;
    
    T(j+1) = a + h*j; %updates by taking another step

    HEFFs1(:,j+1) = -J*s2j;
    HEFFs2(:,j+1) = -J*s1j;      
    
    sDot = dot(s1j,s2j);
    
    E(:,j+1) = -J*sDot; %Energy Hamiltonian
end


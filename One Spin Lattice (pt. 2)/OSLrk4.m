function[T,S1,S2] = OSLrk4(f,g,a,b,s1a,s2a,m,J, HEFFs1, HEFFs2)



h = (b-a)/m; %step size


T = zeros(1,m+1); %time vector (stores all t values)...(1xm+1) matrix

S1 = zeros(3,m+1); %Spin1 vector (stores all Spin1 values)

S2 = zeros(3,m+1); %Spin2 vector 


T(1) = a; %initial time

S1(:,1) = s1a; %intial spin1 (particle1) value
S2(:,1) = s2a; 


for j=1:m, %iterates one time step (entails 4 approximations)

           %like for loop in Java (for(j = 1; j = m; j++))

    

    tj = T(j); %increase j value per iteration of loop

    s1j = S1(:,j);

    s2j = S2(:,j);
    
    heff1 = HEFFs1(:,j);

    heff2 = HEFFs2(:,j);
    

    k1 = h*feval(f,tj,s1j,s2j); %h(run)*(slope of spin1 (s1) tangent line) = approx. rise from last iteration

    l1 = h*feval(g,tj,s1j,s2j); %h(run)*(slope of spin2 (s2) tangent line) = approx. rise from last iteration

    k2 = h*feval(f,tj+(h/2),s1j+(k1/2), s2j+(l1/2)); % "" but using updated values of s1 and 2 (half steps)

    l2 = h*feval(g,tj+(h/2),s1j+(k1/2),s2j+(l1/2)); % "" (half steps)

    k3 = h*feval(f,tj+(h/2),s1j+(k2/2),s2j+(l2/2)); % "" (half steps)

    l3 = h*feval(g,tj+(h/2),s1j+(k2/2),s2j+(l2/2)); % "" (half steps)

    k4 = h*feval(f,tj+h,s1j+k3,s2j+l3); % ""

    l4 = h*feval(g,tj+h,s1j+k3,s2j+l3); % ""

    
    %takes old value of S1 and adds new value to it to update it
    S1(:,j+1) = s1j + (k1 + 2*k2 + 2*k3 + k4)/6; %computes new y value by averaging all 4 approximations of y

    S2(:,j+1) = s2j + (l1 + 2*l2 + 2*l3 + l4)/6; %"" z
    
    T(j+1) = a + h*j; %updates by taking another step

    HEFFs1(:,j+1) = -J*S2(:,j+1); %initial effective h of particle 1
    HEFFs2(:,j+1) = -J*S1(:,j+1);
end


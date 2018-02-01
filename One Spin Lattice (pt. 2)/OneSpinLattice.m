%spins code

%modeling the change in the first spin based SOLELY on the z-spin
%of second atom

%defining variables
spin = zeros (3,2); %holding values for both spins
J = 1; %interaction parameter 

m = 1500; %num steps
a = 0; %intial time
b = 50; %terminal time (seconds)

%to fill the matrix elements of spin1 (random num between 0 and 1)
for i = 1:3, 
        spin(i,1)=randi([0,1]);
end

%to fill x/y-component of spin 2
for i=1:2
    spin(i,2) = randi([0,1]);
end

%function to normalize the matrix columns, so length of vector adds to 1
normc(spin);

s1a = spin(:,[1]); %setting initial particle1 spin = column 1 of spin matrix
s2a = spin(:,[2]);

HEFFs1 = zeros(3,m+1); %h effective vector
HEFFs2 = zeros(3,m+1);

HEFFs1(:,[1]) = J*s2a; %initial effective h of particle 1
HEFFs2(:,[1]) = J*s1a;

heff1 = HEFFs1(:,[1]);
heff2 = HEFFs2(:,[1]);

x = cross(s1a,heff1);
y = cross(s2a,heff2);

f = @(t,s1a,s2a) x; %spin1 function
g = @(t,s2a,s1a) y;

%dSi/dt = -Si x h(eff)i,j
[T,S1,S2] = OSLrk4(f,g,a,b,s1a,s2a,m,J, HEFFs1, HEFFs2);

%plot(T,S1)
%plot(T,S2)
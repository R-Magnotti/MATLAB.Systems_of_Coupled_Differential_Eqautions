%Code for two-spin lattice 
%Makes call to subroutine TSLrk4

%defining variables
spin = zeros (3,2);
J = 1; %interaction parameter 
%M; %magnetization

m = 10000; %num steps
a = 0; %intial time
b = 100.0; %terminal time (seconds)

%to fill the matrix elements (random num between 0 and 1)
for i = 1:3, 
    for j = 1:2,
        spin(i,j)=rand;
    end
end

%function to normalize the matrix columns, so length of vector adds to 1
spin=normc(spin);

s1a = spin(:,1); %setting initial particle1 spin = column 1 of spin matrix
s2a = spin(:,2);

heff1 = J*s2a;
heff2 = J*s1a;
    
f = @(t,heff1,s1a) cross(heff1,s1a); %spin1 function
g = @(t,heff2,s2a) cross(heff2,s2a);

%dSi/dt = -Si x h(eff)i,j
[T,S1,S2,E,HEFFs1,HEFFs2] = TSLrk4(f,g,a,b,s1a,s2a,m,J, heff1, heff2);

%plot(T,E)
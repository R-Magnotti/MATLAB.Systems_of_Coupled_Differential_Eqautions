%Code for three-spin lattice 
%Makes call to subroutine ThSLrk4

%defining variables
spin = zeros(3,3);
J = 1; %interaction parameter 
%M; %magnetization

m = 5000; %num steps
a = 0; %intial time
b = 100.0; %terminal time (seconds)

%to fill the matrix elements (random num between 0 and 1)
for i = 1:3, 
    for j = 1:3,
        spin(i,j)=rand;
    end
end

%function to normalize the matrix columns, so length of vector adds to 1
spin=normc(spin);

s1a = spin(:,1); %setting initial particle1 spin = column 1 of spin matrix
s2a = spin(:,2);
s3a = spin(:,3);

%inc. number of heff by one 
heff1 = -J*(s2a + s3a);
heff2 = -J*(s1a + s3a);
heff3 = -J*(s1a + s2a);
    
f1 = @(t,heff1,s1a) cross(heff1,s1a); %spin1 function
f2 = @(t,heff2,s2a) cross(heff2,s2a);
f3 = @(t,heff3,s3a) cross(heff3,s3a);

%dSi/dt = -Si x h(eff)i,j
[T,S1,S2,S3,E,HEFFs1,HEFFs2,HEFFs3] = ThSLrk4(f1,f2,f3,a,b,s1a,s2a,s3a,m,J,heff1,heff2,heff3);

%plot(T,E)
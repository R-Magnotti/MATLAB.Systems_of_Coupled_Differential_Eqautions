%Code for all-spin lattice 
%Makes call to subroutine ASLrk4

%this code is designed to accept 'n' amount of spins (as opposed to 2 or 3
%and numerically integrate via Runge-Kutta
%However, we were only able to start it int the time allotted

%defining variables
s = 50; %number of spins

spin = zeros (3,s);
J = 1; %interaction parameter 
%M; %magnetization

spinMat = zeros(3,s); %matrix of spin initial values
heffMat = zeros(3,s); %matrix of feff initial values
fMat = zeros(1,s); %matrix of funcions

m = 5000; %num steps
a = 0; %intial time
b = 50.0; %terminal time (seconds)

%to fill the matrix elements (random num between 0 and 1)
for i = 1:3, 
    for j = 1:s,
        spin(i,j)=rand;
    end
end

%function to normalize the matrix columns, so length of vector adds to 1
spin=normc(spin);

%matrix of initial values
for j=1:s,
   
    spinMat(:,j) = spin(:,j); %setting initial particle1 spin = column 1 of spin matrix
    
end

%matrix of intial heff values
for j=1:s,
    
    heffMat(:,j) = -J*(spinMat(:,j-1) + spinMat(:,j+1));
   
end

for j=1:s,
    
    fMat(1,j) = @(t,heffMat(:,j),spinMat(:,j)) cross(heffMat(:,j),spinMat(:,j)); %spin1 function

end

%left off here -- replace initial values with initial values from matrices
%dSi/dt = -Si x h(eff)i,j
[T,S,E,HEFF] = ASLrk4(fMat,a,b,spinMat,m,J,heffMat,s);

%plot(T,E)
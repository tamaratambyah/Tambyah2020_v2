function conc = tridi(a,b,c,d)
% FUNCTION TRIDI 
%
% solves the tridiagonal matrix system using the Thomas algorithm. Returns
% the solution to the tridiagonal system. 

nodes = size(b,1);

% initialise solution 
conc = zeros(nodes, 1);

% pass diagonal and rhs 
bb = b;
dd = d;

for i = 2:nodes
    ff = a(i)/bb(i-1);
    bb(i) = bb(i) - c(i-1)*ff;
    dd(i) = dd(i) - dd(i-1)*ff;
end

% perform back substitution
conc(nodes) = dd(nodes)/bb(nodes);
for i = 1:nodes-1
    j = nodes-i;
    conc(j) = (dd(j)-c(j)*conc(j+1))/bb(j);
end


end
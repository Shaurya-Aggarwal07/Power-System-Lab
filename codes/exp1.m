clc; clear;

% 1 based indexing
% line data: (from-to; R-X; B/2) 
ld = [
    1 2 0.02 0.06 0.03
    1 3 0.08 0.024 0.025
    2 3 0.06 0.018 0.02
];

nb = 3; 
Y = zeros(nb,nb);
% ----------------------------------
nl = 3;
for k=1:nl
    i = ld(k,1); j = ld(k,2);
    r = ld(k,3); x = ld(k,4);
    b = ld(k,5);

    z = r+1j*x;
    y = 1/z;
    bsh = 1j*b;

    Y(i,j) = Y(i,j) - y;
    Y(j,i) = Y(i,j);

    Y(i,i) = Y(i,i) + y + bsh;
    Y(j,j) = Y(j,j) + y + bsh;
end

disp('Ybus matrix');
disp(Y);
% ------------------------------
elim = 3; % bus to be eliminated
keep = setdiff(1:nb, elim); % gives [1 2]

Yaa = Y(keep,keep);
Yab = Y(keep,elim);
Yba = Y(elim,keep);
Ybb = Y(elim,elim);

Yred = Yaa - Yab*(1/Ybb)*Yba; % important

disp('reduced Ybus is:');
disp(Yred);
% -----------------------------
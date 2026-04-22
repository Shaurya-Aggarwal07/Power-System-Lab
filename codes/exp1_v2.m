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

nl = 3;
for k=1:nl
    i = ld(k,1); j = ld(k,2);
    r = ld(k,3); x = ld(k,4);
    b = ld(k,5);

    z = r+1j*x;
    y = 1/z;
    bsh = 1j*b;
    
    % off diagonal
    Y(i,j) = Y(i,j) - y;
    Y(j,i) = Y(i,j);

    % diagonal 
    Y(i,i) = Y(i,i) + y + bsh;
    Y(j,j) = Y(j,j) + y + bsh;
end

disp('Ybus matrix');
disp(Y);
% ------------------------------------------------------------------------
elim = 3;     % bus to be eliminated
Yred = zeros(nb-1, nb-1);

a=1;
for i= 1:nb
    if i == elim
        continue;
    end
    %-----------------
    b= 1;
    for j = 1:nb
        if j == elim
            continue;
        end
        %-------------
        Yred(a,b) = Y(i,j) - (Y(i,elim) * Y(elim,j)) / Y(elim,elim);
        b= b+1;
    end
    %-----------------
    a= a+1;
end
%-----------
disp('reduced Ybus is:');
disp(Yred);
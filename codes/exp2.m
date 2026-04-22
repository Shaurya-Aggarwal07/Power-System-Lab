Y=zeros(3,3);
y=zeros(3,3);
mat=zeros(3,3);

n=length(mat);
% ---------------------------
% (Y/2) wala: shunt admittances

Y(1,2)=0.01i; 
Y(1,3)=0.01i;
Y(2,3)=0.01i;

Y(2,1)=Y(1,2);
Y(3,1)=Y(1,3);
Y(3,2)=Y(2,3);
% ---------------------------
% X wala: series admittances

y(1,2)= 1/(0.1i);
y(1,3)= 1/(0.1i);
y(2,3)= 1/(0.1i);

y(2,1)= y(1,2);
y(3,1)= y(1,3);
y(3,2)= y(2,3);
% -------------------------
for i=1:n
    for j=1:n
        % --------
        if(i==j)
            for k=1:n
                if(k~=i)
                    mat(i,j)= mat(i,j)+ Y(i,k)+y(i,k);
                end
            end
        % --------
        else
            mat(i,j)= -y(i,j);
        end
    end
end

disp('Ybus matrix:');
disp(mat)
% ------------------------------------------------------------------
N=3;
bus={'Slack','PV','PQ'};

Psch =  [0, 0.6661, -2.8653]; % total given: P = Pg-Pd 
Qsch =  [0, 0, -1.2244];  % total given: Q = Qg-Qd 
Qmin  = [0, 0.2, 0];
Qmax  = [0, 2, 0];

V = ones(N,1);  % 3x1 matrix
Vspec = [1.0, 1.05, 1.0];   
V(1) = Vspec(1);   
V(2) = Vspec(2);   

% V will store new values at any instant and Vprev will store old values

eps = 0.001; % epsilon
alpha = 1.6;

maxIter = 50; 
cnt = NaN; % will store total iterations required to converge
% -----------------------------------------------------------------
for iter = 1:maxIter
    Vprev = V;
    % -------------------------
    for i = 2:N 
        if strcmp(bus{i},'PQ')
            
            sumYV=0;
            for j=1:i-1
                sumYV = sumYV + mat(i,j)*V(j);       
            end
            for j=i+1:N
                sumYV = sumYV + mat(i,j)*Vprev(j);   
            end

            Vgs = ( ( Psch(i)-1j*Qsch(i))/conj(V(i) ) - sumYV ) / mat(i,i);
            V(i) = V(i) + alpha*(Vgs - V(i));
        % --------------------
        elseif strcmp(bus{i},'PV')
            %---------
            sumYV_forQ=0;
            for j=1:i-1
                sumYV_forQ = sumYV_forQ + mat(i,j)*V(j);
            end
            for j=i:N
                sumYV_forQ = sumYV_forQ + mat(i,j)*Vprev(j);
            end
           
            Qi = -imag( conj(V(i)) * (sumYV_forQ) );
            if Qi < Qmin(i)
                Qi = Qmin(i);
            elseif Qi > Qmax(i)
                Qi = Qmax(i); 
            end
            %--------
            sumYV=0;
            for j=1:i-1
                sumYV = sumYV + mat(i,j)*V(j);
            end
            for j=i+1:N
                sumYV = sumYV + mat(i,j)*Vprev(j);
            end
            %-----
            Vgs = ( (Psch(i)-1j*Qi)/conj(V(i)) - sumYV ) / mat(i,i);
            V(i) = V(i) + alpha*(Vgs - V(i));
            %----
            V(i) = Vspec(i) * V(i)/abs(V(i));  
        end
    end
    % -------------------------------------------------------
    if max(abs(V - Vprev)) < eps
        cnt = iter;
        break;
    end
end
% ---------------------------------------------------------------
disp(['Converged in ', num2str(cnt), ' iterations']);

disp('Final Bus Voltages:');
disp(V);
% ------------------------------
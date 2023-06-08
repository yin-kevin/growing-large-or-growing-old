function [cov_mat_p] = GMM_std_errors(t_vector_cap,a_vector_cap,p,W_choice_cap)
% This function Computes GMM standard errors using procedure described in Giuseppe's September 10 2021 email

% Remark: this is only valid for the simplified model in growth rates

% Takes as input t_vector_cap and a_vector_cap 
% (these are two auxiliary objects that define which moments we are using),
% the estimates for the four parameters p=[rho_u,sigma2_e,sigma2_z,psi],
% and the weighting matrix chosen (and cut appropriately)

%returns the estimates of the standard errors for each of the four parameters

% rename variables for convenience
rho_u=p(1);
sigma2_e=p(2);
sigma2_z=p(3);
psi=p(4);

t=t_vector_cap;
a=a_vector_cap;
W=W_choice_cap;

% form a matrix G which has as many rows N as there are moments (the size of the vector m) 
% and as many columns K as there are parameters (the size of the vector p). Each entry of G in row j=1,2...N, column i=1,2...K is the derivative of the j-th moment with respect to the i-th parameter. So populating matrix G takes to calculate the derivative of the variance and covariance expressions with respect to each of the parameters, and to evaluate it at \hat p. 

m_size=size(t_vector_cap,1);
p_size=size(p,2);

% We create a matrix G with as many rows as moments and as many columns as parameters
G=zeros(m_size,p_size);

% We populate matrix G using in position (j,i) the derivatives of the theoretical moment j
% with respect to the parameter i, evaluated at p

% Derivative with respect to rho_u (first column)
G(:,1)=(2*(t-1)+a).*rho_u.^(2*(t-1)+a-1).*(1-rho_u).^2*psi;
G(:,1)=G(:,1)+2*rho_u.^(2*(t-1)+a).*(1-rho_u)*psi;

% Derivative with respect to sigma2_e (second column)
G(:,2)=1*(a==0); 
% Remark: the 1* is needed to prevent the creation of a logical array

% Derivative with respect to sigma2_z (third column)
G(:,3)=2*(a==0)-1*(a==1);

% Derivative with respect to psi (fourth column)
G(:,4)=rho_u.^(2*(t-1)+a).*(1-rho_u).^2;

% We have now populated G

% Asymptotic covariance matrix of p is (G'W^{-1} G)^{-1}
cov_mat_p=(G'*(W^(-1))*G)^(-1);
end


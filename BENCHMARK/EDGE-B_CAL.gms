* Calibration of EDGE-BENCHMARK 

* Calibration allows for differences in separation rate and recruitment costs
* Recruiter ratio fixed at 24.99
* Abatement is assumed to be zero in the no-policy steady state


parameters
pbar_ss
lambda_ss
;

Variables
p_ss(ind)
ntot_ss
n_ss(ind)
u_ss(ind)
ubar_ss
theta_ss(ind,ind2)
thetaind_ss(ind)
thetabar_ss
q_ss(ind)
mu(ind)
phi_ss(ind,ind2)
v_n_ss(ind)
v_u_ss(ind)
u_n_ss(ind)
u_u_ss(ind)

y_ss(ind)
vbar_ss(ind)
j_n_ss(ind)
w_ss(ind)
h_ss(ind)
A(ind)
psi
b0
;

Equations
eqn_ntot_ss
eqn_ubar_ss
eqn_u_ss(ind)
eqn_theta_ss(ind,ind2)
eqn_thetaind_ss(ind)
eqn_thetabar_ss
eqn_mu(ind)
eqn_phi_ss(ind,ind2)
eqn_v_n_ss(ind)
eqn_v_u_ss(ind)
eqn1_u_n_ss(ind)
eqn2_u_n_ss(ind)
eqn1_u_u_ss(ind)
eqn2_u_u_ss(ind)
eqn_y_ss(ind)
eqn_vbar_ss(ind)
eqn_j_n_ss(ind)
eqn_w_ss(ind)
eqn_nash_w_ss(ind)
eqn_nash_h_ss(ind)
eqn_mc_ss(ind)
;

* Define fixed Parameters

pbar_ss = 1;
lambda_ss = (c_ss**(-sigma))/pbar_ss;



* Labor Market Equations
eqn_ntot_ss.. ntot_ss =e= sum(ind,n_ss(ind));
eqn_ubar_ss.. ubar_ss =e= sum(ind,u_ss(ind));
eqn_u_ss(ind).. sum(ind2,phi_ss(ind,ind2))*u_ss(ind) =e= sum(ind2,u_ss(ind2)*phi_ss(ind2,ind));

eqn_theta_ss(ind,ind2).. theta_ss(ind,ind2) =e= (vbar_ss(ind)*n_ss(ind)*h_ss(ind))/u_ss(ind2);
eqn_thetaind_ss(ind).. thetaind_ss(ind) =e= (vbar_ss(ind)*n_ss(ind)*h_ss(ind))/sum(ind2,u_ss(ind2));
eqn_thetabar_ss.. thetabar_ss =e= sum(ind,thetaind_ss(ind));
eqn_mu(ind).. mu(ind) =e= q_ss(ind)/(xi*thetabar_ss**(-gamma)+(1-xi)*theta_ss(ind,ind)**(-gamma));
eqn_phi_ss(ind,ind2).. phi_ss(ind,ind2) =e= mu(ind2)*(xi*thetaind_ss(ind2)*thetabar_ss**(-gamma) + (1-xi)*delta(ind,ind2)*theta_ss(ind,ind2)**(1-gamma));
*eqn_phi_ss(ind)..  pi(ind)*n_ss(ind) =e= 

* Household Equations

eqn1_u_n_ss(ind).. u_n_ss(ind) =e= log(c_ss) - (psi/(1+chi))*h_ss(ind)**(1+chi);
eqn2_u_n_ss(ind).. u_n_ss(ind) =e= (1/(1-sigma))*c_ss**(1-sigma) - (psi/(1+chi))*h_ss(ind)**(1+chi);
eqn1_u_u_ss(ind).. u_u_ss(ind) =e= log(c_ss);
eqn2_u_u_ss(ind).. u_u_ss(ind) =e= (1/(1-sigma))*c_ss**(1-sigma);

eqn_v_n_ss(ind).. v_n_ss(ind) =e= ( u_n_ss(ind) + lambda_ss*(1-tau_l0)*w_ss(ind)*h_ss(ind) + pi(ind)*beta*v_u_ss(ind))/(1-beta*(1-pi(ind)));
eqn_v_u_ss(ind).. v_u_ss(ind) =e= ( u_u_ss(ind) + lambda_ss*b0*pbar_ss + beta*sum(ind2,phi_ss(ind,ind2)*v_n_ss(ind2)))/(1-beta*(1-sum(ind2,phi_ss(ind,ind2))));


* Firm Equations

eqn_y_ss(ind)..  y_ss(ind) =e= A(ind)*n_ss(ind)*(1-vbar_ss(ind))*h_ss(ind);
eqn_vbar_ss(ind).. vbar_ss(ind) =e= pi(ind)/(q_ss(ind)*h_ss(ind));
eqn_j_n_ss(ind).. j_n_ss(ind) =e= (p_ss(ind)*A(ind))/(beta*q_ss(ind));
eqn_w_ss(ind).. w_ss(ind) =e= (p_ss(ind)*A(ind)*h_ss(ind) - j_n_ss(ind)*(1-beta*(1-pi(ind))))/(h_ss(ind)*(1+tau_p0));


* Nash Bargaining Equations

eqn_nash_w_ss(ind).. j_n_ss(ind) =e= (eta/(1-eta))*((1+tau_p0)/(1-tau_l0))*(1/lambda_ss)*(v_n_ss(ind)-v_u_ss(ind));
eqn_nash_h_ss(ind).. h_ss(ind) =e= ((((1-tau_l0)/(1+tau_p0))*lambda_ss*A(ind)*p_ss(ind))/psi)**(1/chi);

* Market Clearing

eqn_mc_ss(ind).. y_ss(ind) =e= cd_ss(ind);


* Initialization and fix variables

p_ss.fx(ind) = 1;
p_ss.fx('c') = 1;

* Recruiter productivity fixed 

q_ss.fx('c') = 24.99;
q_ss.fx('d') = 24.99;


* Baseline unemployment set to 5 percent

ntot_ss.fx = 0.95;
ubar_ss.l = 1-ntot_ss.l;
n_ss.l(ind) = ntot_ss.l*cd_ss(ind)/sum(ind2,cd_ss(ind2));
u_ss.l(ind) = (1-ntot_ss.l)*cd_ss(ind)/sum(ind2,cd_ss(ind2));

* Clean hours fixed to 1/3 of time

h_ss.l(ind) = 1/3;
h_ss.fx('c') = 1/3;


* Guess values of employment, disutility of work, unemployment benefits,

psi.l = 6.13;
b0.l = .2758;

v_n_ss.l(ind) = 86.91;
v_u_ss.l(ind) = 86.825;
u_n_ss.l(ind)$(sigma = 1) = log(c_ss) - (psi.l/(1+chi))*h_ss.l(ind)**(1+chi);
u_n_ss.l(ind)$(sigma ne 1) = (1/(1-sigma))*c_ss**(1-sigma) - (psi.l/(1+chi))*h_ss.l(ind)**(1+chi);
u_u_ss.l(ind)$(sigma = 1) = log(c_ss) ;
u_u_ss.l(ind)$(sigma ne 1) = (1/(1-sigma))*c_ss**(1-sigma) ;




* Use model equations to derive initialized values of remaining variables.


vbar_ss.l(ind) = pi(ind)/(q_ss.l(ind)*h_ss.l(ind));
theta_ss.l(ind,ind2) = (vbar_ss.l(ind)*n_ss.l(ind)*h_ss.l(ind))/( u_ss.l(ind2) );
thetaind_ss.l(ind) = sum(ind2,theta_ss.l(ind,ind2));
thetabar_ss.l = sum(ind,thetaind_ss.l(ind));
mu.l(ind) = q_ss.l(ind)/(xi*thetabar_ss.l**(-gamma)+(1-xi)*theta_ss.l(ind,ind)**(-gamma));
phi_ss.l(ind,ind2) = mu.l(ind2)*(xi*thetaind_ss.l(ind2)*thetabar_ss.l**(-gamma) + (1-xi)*delta(ind,ind2)*theta_ss.l(ind,ind2)**(1-gamma));


w_ss.l(ind) =  (v_n_ss.l(ind)*(1-beta*(1-pi(ind))) - u_n_ss.l(ind) - pi(ind)*v_u_ss.l(ind))/(lambda_ss*(1-tau_l0)*h_ss.l(ind));
j_n_ss.l(ind) = (eta/(1-eta))*((1+tau_p0)/(1-tau_l0))*(1/lambda_ss)*(v_n_ss.l(ind)-v_u_ss.l(ind));
y_ss.l(ind) = cd_ss(ind); 

A.l(ind) = y_ss.l(ind)/(n_ss.l(ind)*(1-vbar_ss.l(ind))*h_ss.l(ind));




model edge_b1_cal /
eqn_ntot_ss
eqn_ubar_ss
eqn_u_ss
eqn_theta_ss
eqn_thetaind_ss
eqn_thetabar_ss
eqn_mu
eqn_phi_ss
eqn_v_n_ss
eqn_v_u_ss
eqn1_u_n_ss
eqn1_u_u_ss
eqn_y_ss
eqn_vbar_ss
eqn_j_n_ss
eqn_w_ss
eqn_nash_w_ss
eqn_nash_h_ss
eqn_mc_ss
/
;

model edge_b2_cal /
eqn_ntot_ss
eqn_ubar_ss
eqn_u_ss
eqn_theta_ss
eqn_thetaind_ss
eqn_thetabar_ss
eqn_mu
eqn_phi_ss
eqn_v_n_ss
eqn_v_u_ss
eqn2_u_n_ss
eqn2_u_u_ss
eqn_y_ss
eqn_vbar_ss
eqn_j_n_ss
eqn_w_ss
eqn_nash_w_ss
eqn_nash_h_ss
eqn_mc_ss
/
;

if(sigma = 1,
solve edge_b1_cal using mcp;
);

if(sigma ne 1,
solve edge_b2_cal using mcp;
);

parameter
tlump_ss
prof_ss(ind)
bc_ss
phibar_ss(ind)
pay_ss
w_ss_db(ind)
v_nu_ss
pvlump_ss
;



tlump_ss = (tau_l0+tau_p0)*(sum(ind,n_ss.l(ind)*w_ss.l(ind)*h_ss.l(ind))) - (1-ntot_ss.l)*pbar_ss*b0.l;
prof_ss(ind) = p_ss.l(ind)*y_ss.l(ind) - (1+tau_p0)*w_ss.l(ind)*n_ss.l(ind)*h_ss.l(ind);
bc_ss = (1-tau_l0)*sum(ind,n_ss.l(ind)*w_ss.l(ind)*h_ss.l(ind))  + (1-ntot_ss.l)*pbar_ss*b0.l + tlump_ss + sum(ind,prof_ss(ind)) - pbar_ss*c_ss;
phibar_ss(ind) = beta*sum(ind2,phi_ss.l(ind,ind2)*(v_n_ss.l(ind2) - v_u_ss.l(ind)));

pay_ss(ind) = ((1-tau_l0)/(1+tau_p0))*(1-eta)*( A.l(ind)*p_ss.l(ind)*h_ss.l(ind))
 + eta*((psi.l*(h_ss.l(ind)**(1+chi))/(1+chi))/lambda_ss + b0.l*pbar_ss + phibar_ss(ind)/lambda_ss);

w_ss_db(ind) = pay_ss(ind)/((1-tau_l0)*h_ss.l(ind));

v_nu_ss(ind) = v_n_ss.l(ind) - v_u_ss.l(ind);

pvlump_ss = tlump_ss/(1-beta);


display
prof_ss
bc_ss
psi.l
b0.l
w_ss.l
w_ss_db
vbar_ss.l
v_n_ss.l
v_u_ss.l
v_nu_ss
theta_ss.l
mu.l
phi_ss.l
q_ss.l
ubar_ss.l
phibar_ss
n_ss.l

;


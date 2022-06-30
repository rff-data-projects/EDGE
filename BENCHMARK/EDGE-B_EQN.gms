* Employment Dynamics in General Equilibrium Benchmark Model - EDGE-BENCHMARK

* Two Sector Model with Search Frictions
* U(c,h) = c^(1-sigma)/(1-sigma) - psi h**(1+chi)/(1+chi);
* Numeraire is price of clean good
* Length of period is one month
* No Capital
* Separation and recruitment costs allowed to vary by sector

* Preexisting taxes
* Sectoral Switching Costs and Sticky Wages

* Model code for emissions pricing policies (p_e)


* Define tax variables 
positive variables
tau_l_ss_ct
tau_p_ss_ct
tau_l(t)
tau_p(t)
delta_tau_l
delta_tau_p
;

equations
eqn_tau_l_ss_ct
eqn_tau_p_ss_ct
eqn_tau_l(tp)
eqn_tau_p(tp)
;

eqn_tau_l_ss_ct.. tau_l_ss_ct =e= tau_l0*delta_tau_l;
eqn_tau_p_ss_ct.. tau_p_ss_ct =e= tau_p0*delta_tau_p;
eqn_tau_l(tp).. tau_l(tp) =e= tau_l0*delta_tau_l;
eqn_tau_p(tp).. tau_p(tp) =e= tau_p0*delta_tau_p;

delta_tau_l.l = 1;
delta_tau_p.l = 1;
tau_l_ss_ct.l = tau_l0*delta_tau_l.l;
tau_p_ss_ct.l = tau_p0*delta_tau_p.l;
tau_l.l(tp) = tau_l0*delta_tau_l.l;
tau_p.l(tp) = tau_p0*delta_tau_p.l;
tau_l.fx(t)$(not tp(t)) = tau_l0;
tau_p.fx(t)$(not tp(t)) = tau_p0;

* Steady State

variables
p_ss_ct(ind)
pbar_ss_ct
pnet_ss_ct(ind)

c_ss_ct
lambda_ss_ct
cd_ss_ct(ind)

u_n_ss_ct(ind)
u_u_ss_ct(ind)
v_n_ss_ct(ind)
v_u_ss_ct(ind)
phibar_ss_ct(ind)


w_ss_ct(ind)
h_ss_ct(ind)

vbar_ss_ct(ind)
j_n_ss_ct(ind)
nu_ss_ct(ind)
zbar_ss_ct(ind)
nuprime_ss_ct(ind)
y_ss_ct(ind)
prof_ss_ct(ind)

phi_ss_ct(ind,ind2)
q_ss_ct(ind)
theta_ss_ct(ind,ind2)
thetaind_ss_ct(ind)
thetabar_ss_ct
n_ss_ct(ind)
u_ss_ct(ind)
ntot_ss_ct
ubar_ss_ct

tlump_ss_ct

;


equations
eqn_pbar_ss_ct
eqn_pnet_ss_ct

eqn_hhfoc_ss_ct
eqn_bc_ss_ct
eqn_cd_ss_ct(ind)

eqn1_u_n_ss_ct(ind)
eqn2_u_n_ss_ct(ind)
eqn1_u_u_ss_ct(ind)
eqn2_u_u_ss_ct(ind)
eqn_v_n_ss_ct(ind)
eqn_v_u_ss_ct(ind)
eqn_phibar_ss_ct(ind)

eqn_nash_ss_ct(ind)
eqn_intra_ss_ct(ind)

eqn_firm_ss_ct(ind)
eqn_j_n_ss_ct(ind)
eqn_zbar_ss_ct(ind)
eqn_nu_ss_ct(ind)
eqn_nuprime_ss_ct(ind)
eqn_y_ss_ct(ind)
eqn_prof_ss_ct(ind)

eqn_phi_ss_ct(ind,ind2)
eqn_q_ss_ct(ind)
eqn_theta_ss_ct(ind,ind2)
eqn_thetaind_ss_ct(ind)
eqn_thetabar_ss_ct
eqn_n_ss_ct(ind)
eqn1_u_ss_ct
eqn2_u_ss_ct
eqn_ntot_ss_ct
eqn_ubar_ss_ct

eqn_tlump_ss_ct

eqn_mc_ss_ct(ind)
;

* Aggregate Prices

eqn_pbar_ss_ct.. pbar_ss_ct =e= gammac**(sigmac/(1-sigmac))*(sum(ind,(alphac(ind)**sigmac)*(p_ss_ct(ind))**(1-sigmac)))**(1/(1-sigmac));

eqn_pnet_ss_ct(ind).. pnet_ss_ct(ind) =e= p_ss_ct(ind)-p_e_ss*mu_e(ind)*(1-nu_ss_ct(ind));

* Household

eqn_hhfoc_ss_ct.. c_ss_ct**(-sigma) =e= pbar_ss_ct*lambda_ss_ct;
eqn_bc_ss_ct.. (1-tau_l0)*sum(ind,n_ss_ct(ind)*w_ss_ct(ind)*h_ss_ct(ind)) + ubar_ss_ct*pbar_ss_ct*b0 + pbar_ss_ct*tlump_ss_ct  + sum(ind,prof_ss_ct(ind))  - pbar_ss_ct*c_ss_ct =e= 0;
eqn_cd_ss_ct(ind).. cd_ss_ct(ind) =e= ((gammac*alphac(ind))**sigmac)*((pbar_ss_ct/(p_ss_ct(ind)))**sigmac)*c_ss_ct;

eqn1_u_n_ss_ct(ind).. u_n_ss_ct(ind) =e= log(c_ss_ct) - (psi/(1+chi))*h_ss_ct(ind)**(1+chi);
eqn2_u_n_ss_ct(ind).. u_n_ss_ct(ind) =e= (1/(1-sigma))*c_ss_ct**(1-sigma) - (psi/(1+chi))*h_ss_ct(ind)**(1+chi);
eqn1_u_u_ss_ct(ind).. u_u_ss_ct(ind) =e= log(c_ss_ct);
eqn2_u_u_ss_ct(ind).. u_u_ss_ct(ind) =e= (1/(1-sigma))*c_ss_ct**(1-sigma);

eqn_v_n_ss_ct(ind).. v_n_ss_ct(ind) =e= ( u_n_ss_ct(ind) + lambda_ss_ct*(1-tau_l_ss_ct)*w_ss_ct(ind)*h_ss_ct(ind) + pi(ind)*beta*v_u_ss_ct(ind))/(1-beta*(1-pi(ind)));
eqn_v_u_ss_ct(ind).. v_u_ss_ct(ind) =e= ( u_u_ss_ct(ind) + lambda_ss_ct*b0*pbar_ss_ct + beta*sum(ind2,phi_ss_ct(ind,ind2)*v_n_ss_ct(ind2)))/(1-beta*(1-sum(ind2,phi_ss_ct(ind,ind2))));
eqn_phibar_ss_ct(ind).. phibar_ss_ct(ind) =e= beta*sum(ind2,phi_ss_ct(ind,ind2)*(v_n_ss_ct(ind2) - v_u_ss_ct(ind)));

* Nash Bargaining

eqn_nash_ss_ct(ind).. (1-tau_l_ss_ct)*w_ss_ct(ind)*h_ss_ct(ind) =e= ((1-tau_l_ss_ct)/(1+tau_p_ss_ct))*(1-eta)*( A(ind)*pnet_ss_ct(ind)*h_ss_ct(ind)*(1-zbar_ss_ct(ind))  )
 + eta*((psi*(h_ss_ct(ind)**(1+chi))/(1+chi))/lambda_ss_ct + b0*pbar_ss_ct + phibar_ss_ct(ind)/lambda_ss_ct);
eqn_intra_ss_ct(ind).. (1-tau_l_ss_ct)*lambda_ss_ct*A(ind)*pnet_ss_ct(ind)*(1-zbar_ss_ct(ind)) =e= (1+tau_p_ss_ct)*psi*(h_ss_ct(ind)**chi);

* Firm

eqn_firm_ss_ct(ind).. A(ind)*pnet_ss_ct(ind)*(1-zbar_ss_ct(ind))  =e= beta*q_ss_ct(ind)*j_n_ss_ct(ind);
eqn_j_n_ss_ct(ind).. j_n_ss_ct(ind) =e= (A(ind)*h_ss_ct(ind)*pnet_ss_ct(ind)*(1-zbar_ss_ct(ind)) - (1+tau_p_ss_ct)*w_ss_ct(ind)*h_ss_ct(ind))/(1-beta*(1-pi(ind)));
eqn_nu_ss_ct(ind).. pnet_ss_ct(ind)*nuprime_ss_ct(ind) =e= p_e_ss*mu_e(ind)*(1 -zbar_ss_ct(ind));
eqn_zbar_ss_ct(ind).. zbar_ss_ct(ind) =e= theta1*nu_ss_ct(ind)**theta2;
eqn_nuprime_ss_ct(ind).. nuprime_ss_ct(ind) =e= theta1*theta2*nu_ss_ct(ind)**(theta2-1);
eqn_y_ss_ct(ind).. y_ss_ct(ind) =e= A(ind)*h_ss_ct(ind)*n_ss_ct(ind)*(1-vbar_ss_ct(ind))*(1-zbar_ss_ct(ind));
eqn_prof_ss_ct(ind).. prof_ss_ct(ind) =e= pnet_ss_ct(ind)*y_ss_ct(ind) - (1+tau_p_ss_ct)*w_ss_ct(ind)*n_ss_ct(ind)*h_ss_ct(ind);


* Labor Market

eqn_phi_ss_ct(ind,ind2).. phi_ss_ct(ind,ind2) =e= mu(ind2)*(xi*thetaind_ss_ct(ind2)*thetabar_ss_ct**(-gamma) + (1-xi)*delta(ind,ind2)*theta_ss_ct(ind,ind2)**(1-gamma));
eqn_q_ss_ct(ind).. q_ss_ct(ind) =e= mu(ind)*(xi*thetabar_ss_ct**(-gamma)+(1-xi)*theta_ss_ct(ind,ind)**(-gamma));
eqn_theta_ss_ct(ind,ind2).. theta_ss_ct(ind,ind2) =e= vbar_ss_ct(ind)*h_ss_ct(ind)*n_ss_ct(ind)/(u_ss_ct(ind2));
eqn_thetaind_ss_ct(ind).. thetaind_ss_ct(ind) =e= vbar_ss_ct(ind)*h_ss_ct(ind)*n_ss_ct(ind)/(sum(ind2,u_ss_ct(ind2)));
eqn_thetabar_ss_ct.. thetabar_ss_ct =e= sum(ind,thetaind_ss_ct(ind));
eqn_n_ss_ct(ind).. pi(ind)*n_ss_ct(ind) =e= sum(ind2,u_ss_ct(ind2)*phi_ss_ct(ind2,ind));
*eqn_u_ss_ct(ind).. sum(ind2,phi_ss_ct(ind,ind2))*u_ss_ct(ind) =e= pi(ind)*n_ss_ct(ind);
eqn1_u_ss_ct.. sum(ind2,phi_ss_ct('c',ind2))*u_ss_ct('c') =e= pi('c')*n_ss_ct('c');
eqn2_u_ss_ct.. u_ss_ct('d') =e= 1 - ntot_ss_ct - u_ss_ct('c');
eqn_ntot_ss_ct.. ntot_ss_ct =e= sum(ind,n_ss_ct(ind));
eqn_ubar_ss_ct.. ubar_ss_ct =e= sum(ind,u_ss_ct(ind));

* Government

eqn_tlump_ss_ct.. tlump_ss_ct =e= (sum(ind,(tau_l_ss_ct+tau_p_ss_ct)*n_ss_ct(ind)*w_ss_ct(ind)*h_ss_ct(ind) + p_e_ss*mu_e(ind)*(1-nu_ss_ct(ind))*y_ss_ct(ind)) - (ubar_ss_ct)*pbar_ss_ct*b0)/pbar_ss_ct;

* Market Clearing

eqn_mc_ss_ct(ind).. cd_ss_ct(ind) =e= y_ss_ct(ind) ;


* Initialize with values from non-policy steady state

p_ss_ct.fx('c') = 1;
p_ss_ct.l('d') = 1;
pbar_ss_ct.l = 1;
pnet_ss_ct.l(ind) = p_ss_ct.l(ind);

c_ss_ct.l = c_ss;
lambda_ss_ct.l = lambda_ss;
cd_ss_ct.l(ind) = cd_ss(ind);

u_n_ss_ct.l(ind) = u_n_ss.l(ind);
u_u_ss_ct.l(ind) = u_u_ss.l(ind);
v_n_ss_ct.l(ind) = v_n_ss.l(ind);
v_u_ss_ct.l(ind) = v_u_ss.l(ind);
phibar_ss_ct.l(ind) = beta*sum(ind2,phi_ss.l(ind,ind2)*(v_n_ss.l(ind2)-v_u_ss.l(ind)));

w_ss_ct.l(ind) = w_ss.l(ind);
h_ss_ct.l(ind) = h_ss.l(ind);

vbar_ss_ct.l(ind) = vbar_ss.l(ind);
j_n_ss_ct.l(ind) = j_n_ss.l(ind);
nu_ss_ct.l('d') = 0;
nu_ss_ct.fx('c') = 0;
zbar_ss_ct.l('d') = 0;
zbar_ss_ct.fx('c') = 0;
nuprime_ss_ct.l(ind) = 0;
y_ss_ct.l(ind) = y_ss.l(ind);
prof_ss_ct.l(ind) = prof_ss(ind);

phi_ss_ct.l(ind,ind2) = phi_ss.l(ind,ind2);
q_ss_ct.l(ind) = q_ss.l(ind);
theta_ss_ct.l(ind,ind2) = theta_ss.l(ind,ind2);
thetaind_ss_ct.l(ind) = thetaind_ss.l(ind);
thetabar_ss_ct.l = thetabar_ss.l;
u_ss_ct.l(ind) = u_ss.l(ind);
n_ss_ct.l(ind) = n_ss.l(ind);
ntot_ss_ct.l = ntot_ss.l;
ubar_ss_ct.l = ubar_ss.l;

tlump_ss_ct.l = tlump_ss;

mu.fx(ind) = mu.l(ind);
A.fx(ind) = A.l(ind);
psi.fx = psi.l;
b0.fx = b0.l;



model edge_b1_ss_ct /
eqn_pbar_ss_ct.pbar_ss_ct
eqn_pnet_ss_ct.pnet_ss_ct
eqn_hhfoc_ss_ct.c_ss_ct
eqn_bc_ss_ct.lambda_ss_ct
eqn_cd_ss_ct.cd_ss_ct
eqn_v_n_ss_ct.v_n_ss_ct
eqn_v_u_ss_ct.v_u_ss_ct
eqn1_u_n_ss_ct.u_n_ss_ct
eqn1_u_u_ss_ct.u_u_ss_ct
eqn_phibar_ss_ct.phibar_ss_ct
eqn_nash_ss_ct.w_ss_ct
eqn_intra_ss_ct.h_ss_ct
eqn_firm_ss_ct.vbar_ss_ct
eqn_j_n_ss_ct.j_n_ss_ct
eqn_nu_ss_ct.nu_ss_ct
eqn_zbar_ss_ct.zbar_ss_ct
eqn_nuprime_ss_ct.nuprime_ss_ct
eqn_y_ss_ct.y_ss_ct
eqn_prof_ss_ct.prof_ss_ct
eqn_phi_ss_ct.phi_ss_ct
eqn_q_ss_ct.q_ss_ct
eqn_theta_ss_ct.theta_ss_ct
eqn_thetaind_ss_ct.thetaind_ss_ct
eqn_thetabar_ss_ct.thetabar_ss_ct
eqn_n_ss_ct.n_ss_ct
eqn1_u_ss_ct
eqn2_u_ss_ct
eqn_ntot_ss_ct.ntot_ss_ct
eqn_ubar_ss_ct.ubar_ss_ct
eqn_tlump_ss_ct.tlump_ss_ct
eqn_mc_ss_ct.p_ss_ct
eqn_tau_l_ss_ct.tau_l_ss_ct
eqn_tau_p_ss_Ct.tau_p_ss_ct

/
;


model edge_b2_ss_ct /
eqn_pbar_ss_ct.pbar_ss_ct
eqn_pnet_ss_ct.pnet_ss_ct
eqn_hhfoc_ss_ct.c_ss_ct
eqn_bc_ss_ct.lambda_ss_ct
eqn_cd_ss_ct.cd_ss_ct
eqn_v_n_ss_ct.v_n_ss_ct
eqn_v_u_ss_ct.v_u_ss_ct
eqn2_u_n_ss_ct.u_n_ss_ct
eqn2_u_u_ss_ct.u_u_ss_ct
eqn_phibar_ss_ct.phibar_ss_ct
eqn_nash_ss_ct.w_ss_ct
eqn_intra_ss_ct.h_ss_ct
eqn_firm_ss_ct.vbar_ss_ct
eqn_j_n_ss_ct.j_n_ss_ct
eqn_nu_ss_ct.nu_ss_ct
eqn_zbar_ss_ct.zbar_ss_ct
eqn_nuprime_ss_ct.nuprime_ss_ct
eqn_y_ss_ct.y_ss_ct
eqn_prof_ss_ct.prof_ss_ct
eqn_phi_ss_ct.phi_ss_ct
eqn_q_ss_ct.q_ss_ct
eqn_theta_ss_ct.theta_ss_ct
eqn_thetaind_ss_ct.thetaind_ss_ct
eqn_thetabar_ss_ct.thetabar_ss_ct
eqn_n_ss_ct.n_ss_ct
eqn1_u_ss_ct
eqn2_u_ss_ct
eqn_ntot_ss_ct.ntot_ss_ct
eqn_ubar_ss_ct.ubar_ss_ct
eqn_tlump_ss_ct.tlump_ss_ct
eqn_mc_ss_ct.p_ss_ct
eqn_tau_l_ss_ct.tau_l_ss_ct
eqn_tau_p_ss_Ct.tau_p_ss_ct

/
;



* Transition Model

parameters
n0(ind)
u0(ind)
lambda_tt
j_n_tt(ind)
v_n_tt(ind)
v_u_tt(ind)
wdelta_tt(ind)
wstar_tt(ind)
p_e(t)
;

variables
pbar(t)
pnet(ind,t)
c(t)
m(t)
cd(ind,t)
u_n(ind,t)
u_u(ind,t)
v_n(ind,t)
v_u(ind,t)
phibar(ind,t)
w0(ind,t)
wdelta(ind,t)
wstar(ind,t)
w(ind,t)
h(ind,t)
vbar(ind,t)
j_n(ind,t)
zbar(ind,t)
nu(ind,t)
nuprime(ind,t)
y(ind,t)
prof(ind,t)
phi(ind,ind2,t)
q(ind,t)
theta(ind,ind2,t)
thetaind(ind,t)
thetabar(t)
n(ind,t)
u(ind,t)
ntot(t)
ubar(t)
tlump(t)
;

positive variables
p(ind,t)
lambda(t)
;


equations
eqn_pbar(t)
eqn_pnet(ind,t)

eqn_hhfoc(t)
eqn_hheuler(t)
eqn_hheuler_tc(t)
eqn_hhbc(t)
eqn_cd(ind,t)

eqn1_u_n(ind,t)
eqn2_u_n(ind,t)
eqn1_u_u(ind,t)
eqn2_u_u(ind,t)
eqn_v_n(ind,t)
eqn_v_n_tc(ind,t)
eqn_v_u(ind,t)
eqn_v_u_tc(ind,t)
eqn_phibar(ind,t)
eqn_phibar_tc(ind,t)

eqn_nash(ind,t)
eqn_wdelta(ind,t)
eqn_wdelta_tc(ind,t)
eqn_wstar(ind,t)
eqn_wstar_tc(ind,t)
eqn_w(ind,t)
eqn_intra(ind,t)

eqn_firmfoc(ind,t)
eqn_firmfoc_tc(ind,t)
eqn_vbarlast(ind,t)
eqn_j_n(ind,t)
eqn_j_n_tc(ind,t)
eqn_nu(ind,t)
eqn_zbar(ind,t)
eqn_nuprime(ind,t)
eqn_y(ind,t)
eqn_prof(ind,t)

eqn_phi(ind,ind2,t)
eqn_q(ind,t)
eqn_theta(ind,ind2,t)
eqn_thetaind(ind,t)
eqn_thetabar(t)
eqn_n(ind,t)
eqn1_u(t)
eqn2_u(t)
eqn_ntot(t)
eqn_ubar

eqn_tlump(t)

eqn_mc(ind,t)
;

* Aggregate Prices

eqn_pbar(t).. pbar(t) =e= gammac**(sigmac/(1-sigmac))*(sum(ind,(alphac(ind)**sigmac)*(p(ind,t))**(1-sigmac)))**(1/(1-sigmac));
eqn_pnet(ind,t).. pnet(ind,t) =e= p(ind,t) - p_e(t)*mu_e(ind)*(1-nu(ind,t)); ;

* Household

eqn_hhfoc(t).. c(t)**(-sigma) =e= pbar(t)*lambda(t);
eqn_hheuler(t)$(not tlast(t)).. m(t) =e= beta*lambda(t+1)/lambda(t);
eqn_hheuler_tc(tlast).. m(tlast) =e= beta*lambda_ss_ct/lambda(tlast);
eqn_hhbc(t).. (1-tau_l(t))*sum(ind,w(ind,t)*n(ind,t)*h(ind,t)) + ubar(t)*pbar(t)*b0 + tlump(t)*pbar(t) + sum(ind,prof(ind,t)) - pbar(t)*c(t) =g= 0;
eqn_cd(ind,t).. cd(ind,t) =e= ((gammac*alphac(ind))**sigmac)*((pbar(t)/(p(ind,t)))**sigmac)*c(t);

eqn1_u_n(ind,t).. u_n(ind,t) =e= log(c(t)) - (psi/(1+chi))*h(ind,t)**(1+chi);
eqn2_u_n(ind,t).. u_n(ind,t) =e= (1/(1-sigma))*c(t)**(1-sigma) - (psi/(1+chi))*h(ind,t)**(1+chi);
eqn1_u_u(ind,t).. u_u(ind,t) =e= log(c(t));
eqn2_u_u(ind,t).. u_u(ind,t) =e= (1/(1-sigma))*c(t)**(1-sigma);

eqn_v_n(ind,t)$(not tlast(t)).. v_n(ind,t) =e= u_n(ind,t) + lambda(t)*(1-tau_l(t))*w(ind,t)*h(ind,t) + pi(ind)*beta*v_u(ind,t+1) + (1-pi(ind))*beta*v_n(ind,t+1) ;
eqn_v_n_tc(ind,tlast).. v_n(ind,tlast) =e= u_n(ind,tlast) + lambda(tlast)*(1-tau_l(tlast))*w(ind,tlast)*h(ind,tlast) + pi(ind)*beta*v_u_ss_ct(ind) + (1-pi(ind))*beta*v_n_ss_ct(ind) ;
eqn_v_u(ind,t)$(not tlast(t)).. v_u(ind,t) =e=  u_u(ind,t) + lambda(t)*b0*pbar(t) + beta*sum(ind2,phi(ind,ind2,t)*v_n(ind2,t+1)) + beta*(1-sum(ind2,phi(ind,ind2,t)))*v_u(ind,t+1);
eqn_v_u_tc(ind,tlast).. v_u(ind,tlast) =e=  u_u(ind,tlast) + lambda(tlast)*b0*pbar(tlast) + beta*sum(ind2,phi(ind,ind2,tlast)*v_n_ss_ct(ind2)) + beta*(1-sum(ind2,phi(ind,ind2,tlast)))*v_u_ss_ct(ind);
eqn_phibar(ind,t)$(not tlast(t)).. phibar(ind,t) =e= beta*sum(ind2,phi(ind,ind2,t)*(v_n(ind2,t+1) - v_u(ind,t+1)));
eqn_phibar_tc(ind,tlast).. phibar(ind,tlast) =e= beta*sum(ind2,phi(ind,ind2,tlast)*(v_n_ss_ct(ind2) - v_u_ss_ct(ind)));


* Nash Bargaining

eqn_nash(ind,t).. w0(ind,t) =e= (1-eta)*((1-tau_l(t))/(1+tau_p(t)))*(A(ind)*pnet(ind,t)*h(ind,t)*(1-zbar(ind,t)) )
 + eta*((psi*(h(ind,t)**(1+chi))/(1+chi))/lambda(t) + b0*pbar(t) + phibar(ind,t)/lambda(t));
eqn_wdelta(ind,t)$(not tlast(t)).. wdelta(ind,t) =e= h(ind,t) + rho*(1-pi(ind))*m(t)*wdelta(ind,t+1);
eqn_wdelta_tc(ind,tlast).. wdelta(ind,tlast) =e= h(ind,tlast) + rho*(1-pi(ind))*m(tlast)*(h_ss_ct(ind)/(1-rho*beta*(1-pi(ind))));
eqn_wstar(ind,t)$(not tlast(t)).. wdelta(ind,t)*(1-tau_l(t))*wstar(ind,t) =e= w0(ind,t) + rho*(1-pi(ind))*m(t)*wdelta(ind,t+1)*(1-tau_l(t+1))*wstar(ind,t+1);
eqn_wstar_tc(ind,tlast).. wdelta(ind,tlast)*(1-tau_l(tlast))*wstar(ind,tlast) =e= w0(ind,tlast) + rho*(1-pi(ind))*m(tlast)*(h_ss_ct(ind)/(1-rho*beta*(1-pi(ind))))*(1-tau_l_ss_ct)*w_ss_ct(ind);
eqn_w(ind,t).. w(ind,t) =e= (1-rho)*wstar(ind,t) + rho*w(ind,t-1) + (rho*w_ss.l(ind))$tfirst(t);

eqn_intra(ind,t).. (1-tau_l(t))*lambda(t)*A(ind)*pnet(ind,t)*(1-zbar(ind,t)) =e= (1+tau_p(t))*psi*(h(ind,t)**chi);

* Firm

eqn_firmfoc(ind,t)$(not tlast(t)).. A(ind)*pnet(ind,t)*(1-zbar(ind,t)) =e= m(t)*q(ind,t)*j_n(ind,t+1);
eqn_firmfoc_tc(ind,tlast).. A(ind)*pnet(ind,tlast)*(1-zbar(ind,tlast)) =e= m(tlast)*q(ind,tlast)*j_n_ss_ct(ind);
eqn_vbarlast(ind,tlast).. vbar(ind,tlast) =e= vbar_ss_ct(ind);
eqn_j_n(ind,t)$(not tlast(t)).. j_n(ind,t) =e= A(ind)*h(ind,t)*pnet(ind,t)*(1-zbar(ind,t)) - (1+tau_p(t))*w(ind,t)*h(ind,t) + (1-pi(ind))*m(t)*j_n(ind,t+1);
eqn_j_n_tc(ind,tlast).. j_n(ind,tlast) =e= A(ind)*h(ind,tlast)*pnet(ind,tlast)*(1-zbar(ind,tlast)) - (1+tau_p(tlast))*w(ind,tlast)*h(ind,tlast) + (1-pi(ind))*m(tlast)*j_n_ss_ct(ind);
eqn_nu(ind,t).. pnet(ind,t)*nuprime(ind,t) =e= p_e(t)*mu_e(ind)*(1-zbar(ind,t));
eqn_zbar(ind,t).. zbar(ind,t) =e= theta1*nu(ind,t)**theta2;
eqn_nuprime(ind,t).. nuprime(ind,t) =e= theta1*theta2*nu(ind,t)**(theta2-1);
eqn_y(ind,t).. y(ind,t) =e= A(ind)*n(ind,t)*h(ind,t)*(1-vbar(ind,t))*(1-zbar(ind,t));
eqn_prof(ind,t).. prof(ind,t) =e= pnet(ind,t)*y(ind,t) - (1+tau_p(t))*w(ind,t)*n(ind,t)*h(ind,t);

* Labor Market

eqn_phi(ind,ind2,t).. phi(ind,ind2,t) =e= mu(ind2)*(xi*thetaind(ind2,t)*thetabar(t)**(-gamma) + (1-xi)*delta(ind,ind2)*theta(ind,ind2,t)**(1-gamma));
eqn_q(ind,t).. q(ind,t) =e= mu(ind)*(xi*thetabar(t)**(-gamma)+(1-xi)*theta(ind,ind,t)**(-gamma));
eqn_theta(ind,ind2,t).. theta(ind,ind2,t) =e= vbar(ind,t)*h(ind,t)*n(ind,t)/(u(ind2,t));
eqn_thetaind(ind,t).. thetaind(ind,t) =e= vbar(ind,t)*h(ind,t)*n(ind,t)/(sum(ind2,u(ind2,t)));
eqn_thetabar(t).. thetabar(t) =e= sum(ind,thetaind(ind,t));
eqn_n(ind,t).. n(ind,t) =e= n_ss.l(ind)$tfirst(t) + (1-pi(ind))*n(ind,t-1) + q(ind,t-1)*vbar(ind,t-1)*h(ind,t-1)*n(ind,t-1);
eqn_ntot(t).. ntot(t) =e= sum(ind,n(ind,t));
eqn1_u(t).. u('c',t) =e= u_ss.l('c')$tfirst(t) + (1-sum(ind2,phi('c',ind2,t-1)))*u('c',t-1) + pi('c')*n('c',t-1);
eqn2_u(t).. u('d',t) =e= 1 - ntot(t) - u('c',t);
eqn_ubar(t).. ubar(t) =e= sum(ind,u(ind,t));

* Government

eqn_tlump(t).. tlump(t) =e= (sum(ind,(tau_l(t)+tau_p(t))*n(ind,t)*w(ind,t)*h(ind,t) + p_e(t)*mu_e(ind)*(1-nu(ind,t))*y(ind,t)) - ubar(t)*pbar(t)*b0)/pbar(t);

* Market Clearing

eqn_mc(ind,t).. y(ind,t) - cd(ind,t) =g= 0;

* Intertemporal revenue neutrality

equation
eqn_diff
;


eqn_diff.. sum(t,tlump(t)*betat(t))+(tlump_ss_ct*(beta**nt) )/(1-beta) =e= pvlump_ss;



* Initialize
p.fx('c',t) = 1;
p.l('d',t) = p_ss_ct.l('d');
pbar.l(t) = pbar_ss_ct.l;
pnet.l(ind,t) = pnet_ss_ct.l(ind);
lambda.l(t) = lambda_ss_ct.l;
c.l(t) = c_ss_ct.l;
cd.l(ind,t) = cd_ss_ct.l(ind);
m.l(t) = beta;
u_n.l(ind,t) = u_n_ss_ct.l(ind);
u_u.l(ind,t) = u_u_ss_ct.l(ind);
v_n.l(ind,t) = v_n_ss_ct.l(ind);
v_u.l(ind,t) = v_u_ss_ct.l(ind);
phibar.l(ind,t) = phibar_ss_ct.l(ind);
w0.l(ind,t) = (1-tau_l0)*w_ss_ct.l(ind)*h_ss_ct.l(ind);
wdelta.l(ind,t) = h_ss_ct.l(ind)/(1-rho*beta*(1-pi(ind)));
wstar.l(ind,t) = w_ss_ct.l(ind);
w.l(ind,t) = w_ss_ct.l(ind);
h.l(ind,t) = h_ss_ct.l(ind);
vbar.l(ind,t) = vbar_ss_ct.l(ind);
j_n.l(ind,t) = j_n_ss_ct.l(ind);
nu.l('d',t) = nu_ss_ct.l('d');
nu.fx('c',t) = 0;
zbar.l('d',t) = zbar_ss_ct.l('d');
zbar.fx('c',t) = 0;
nuprime.l('d',t) = nuprime_ss_ct.l('d');
nuprime.l('c',t) = 0;
y.l(ind,t) = y_ss_ct.l(ind);
prof.l(ind,t) = prof_ss_ct.l(ind);
phi.l(ind,ind2,t) = phi_ss_ct.l(ind,ind2);
q.l(ind,t) = q_ss_ct.l(ind);
theta.l(ind,ind2,t) = theta_ss_ct.l(ind,ind2);
thetaind.l(ind,t) = thetaind_ss_ct.l(ind);
thetabar.l(t) = thetabar_ss_ct.l;
n.l(ind,t) = n_ss_ct.l(ind);
ntot.l(t) = ntot_ss_ct.l;
u.l(ind,t) = u_ss_ct.l(ind);
ubar.l(t) = ubar_ss_ct.l;
tlump.l(t) = tlump_ss_ct.l;


;


vbar.lo(ind,t)$(not tlast(t)) = 0;

model edge_b1_lump /

eqn_pbar_ss_ct.pbar_ss_ct
eqn_pnet_ss_ct.pnet_ss_ct
eqn_hhfoc_ss_ct.c_ss_ct
eqn_bc_ss_ct.lambda_ss_ct
eqn_cd_ss_ct.cd_ss_ct
eqn_v_n_ss_ct.v_n_ss_ct
eqn_v_u_ss_ct.v_u_ss_ct
eqn1_u_n_ss_ct.u_n_ss_ct
eqn1_u_u_ss_ct.u_u_ss_ct
eqn_phibar_ss_ct.phibar_ss_ct
eqn_nash_ss_ct.w_ss_ct
eqn_intra_ss_ct.h_ss_ct
eqn_firm_ss_ct.vbar_ss_ct
eqn_j_n_ss_ct.j_n_ss_ct
eqn_nu_ss_ct.nu_ss_ct
eqn_zbar_ss_ct.zbar_ss_ct
eqn_nuprime_ss_ct.nuprime_ss_ct
eqn_y_ss_ct.y_ss_ct
eqn_prof_ss_ct.prof_ss_ct
eqn_phi_ss_ct.phi_ss_ct
eqn_q_ss_ct.q_ss_ct
eqn_theta_ss_ct.theta_ss_ct
eqn_thetaind_ss_ct.thetaind_ss_ct
eqn_thetabar_ss_ct.thetabar_ss_ct
eqn_n_ss_ct.n_ss_ct
eqn1_u_ss_ct
eqn2_u_ss_ct
eqn_ntot_ss_ct.ntot_ss_ct
eqn_ubar_ss_ct.ubar_ss_ct
eqn_tlump_ss_ct.tlump_ss_ct
eqn_mc_ss_ct.p_ss_ct


eqn_pbar.pbar
eqn_pnet.pnet
eqn_hhfoc.c
eqn_hheuler
eqn_hheuler_tc
eqn_hhbc.lambda
eqn_cd.cd
eqn1_u_n.u_n
eqn1_u_u.u_u
eqn_v_n.v_n
eqn_v_n_tc
eqn_v_u.v_u
eqn_v_u_tc
eqn_phibar.phibar
eqn_phibar_tc
eqn_nash.w0,
eqn_wdelta.wdelta,
eqn_wdelta_tc,
eqn_wstar.wstar,
eqn_wstar_tc,
eqn_w.w,
eqn_intra.h
eqn_firmfoc.vbar
eqn_vbarlast
eqn_j_n.j_n
eqn_j_n_tc
eqn_nu.nu
eqn_zbar.zbar
eqn_nuprime.nuprime
eqn_y.y
eqn_prof.prof
eqn_phi.phi
eqn_q.q
eqn_theta.theta
eqn_thetaind.thetaind
eqn_thetabar.thetabar
eqn_n.n
eqn_ntot.ntot
eqn1_u
eqn2_u
eqn_ubar.ubar
eqn_tlump.tlump
eqn_mc.p

eqn_tau_l_ss_ct.tau_l_ss_ct
eqn_tau_p_ss_ct.tau_p_ss_ct
eqn_tau_l.tau_l
eqn_tau_p.tau_p
/;

model edge_b2_lump /
eqn_pbar_ss_ct.pbar_ss_ct
eqn_pnet_ss_ct.pnet_ss_ct
eqn_hhfoc_ss_ct.c_ss_ct
eqn_bc_ss_ct.lambda_ss_ct
eqn_cd_ss_ct.cd_ss_ct
eqn_v_n_ss_ct.v_n_ss_ct
eqn_v_u_ss_ct.v_u_ss_ct
eqn2_u_n_ss_ct.u_n_ss_ct
eqn2_u_u_ss_ct.u_u_ss_ct
eqn_phibar_ss_ct.phibar_ss_ct
eqn_nash_ss_ct.w_ss_ct
eqn_intra_ss_ct.h_ss_ct
eqn_firm_ss_ct.vbar_ss_ct
eqn_j_n_ss_ct.j_n_ss_ct
eqn_nu_ss_ct.nu_ss_ct
eqn_zbar_ss_ct.zbar_ss_ct
eqn_nuprime_ss_ct.nuprime_ss_ct
eqn_y_ss_ct.y_ss_ct
eqn_prof_ss_ct.prof_ss_ct
eqn_phi_ss_ct.phi_ss_ct
eqn_q_ss_ct.q_ss_ct
eqn_theta_ss_ct.theta_ss_ct
eqn_thetaind_ss_ct.thetaind_ss_ct
eqn_thetabar_ss_ct.thetabar_ss_ct
eqn_n_ss_ct.n_ss_ct
eqn1_u_ss_ct
eqn2_u_ss_ct
eqn_ntot_ss_ct.ntot_ss_ct
eqn_ubar_ss_ct.ubar_ss_ct
eqn_tlump_ss_ct.tlump_ss_ct
eqn_mc_ss_ct.p_ss_ct

eqn_pbar.pbar
eqn_pnet.pnet
eqn_hhfoc.c
eqn_hheuler
eqn_hheuler_tc
eqn_hhbc.lambda
eqn_cd.cd
eqn2_u_n.u_n
eqn2_u_u.u_u
eqn_v_n.v_n
eqn_v_n_tc
eqn_v_u.v_u
eqn_v_u_tc
eqn_phibar.phibar
eqn_phibar_tc
eqn_nash.w0,
eqn_wdelta.wdelta,
eqn_wdelta_tc,
eqn_wstar.wstar,
eqn_wstar_tc,
eqn_w.w,
eqn_intra.h
eqn_firmfoc.vbar
eqn_vbarlast
eqn_j_n.j_n
eqn_j_n_tc
eqn_nu.nu
eqn_zbar.zbar
eqn_nuprime.nuprime
eqn_y.y
eqn_prof.prof
eqn_phi.phi
eqn_q.q
eqn_theta.theta
eqn_thetaind.thetaind
eqn_thetabar.thetabar
eqn_n.n
eqn_ntot.ntot
eqn1_u
eqn2_u
eqn_ubar.ubar
eqn_tlump.tlump
eqn_mc.p

eqn_tau_l_ss_ct.tau_l_ss_ct
eqn_tau_p_ss_ct.tau_p_ss_ct
eqn_tau_l.tau_l
eqn_tau_p.tau_p
/;

model edge_b1_taul /
eqn_pbar_ss_ct.pbar_ss_ct
eqn_pnet_ss_ct.pnet_ss_ct
eqn_hhfoc_ss_ct.c_ss_ct
eqn_bc_ss_ct.lambda_ss_ct
eqn_cd_ss_ct.cd_ss_ct
eqn_v_n_ss_ct.v_n_ss_ct
eqn_v_u_ss_ct.v_u_ss_ct
eqn1_u_n_ss_ct.u_n_ss_ct
eqn1_u_u_ss_ct.u_u_ss_ct
eqn_phibar_ss_ct.phibar_ss_ct
eqn_nash_ss_ct.w_ss_ct
eqn_intra_ss_ct.h_ss_ct
eqn_firm_ss_ct.vbar_ss_ct
eqn_j_n_ss_ct.j_n_ss_ct
eqn_nu_ss_ct.nu_ss_ct
eqn_zbar_ss_ct.zbar_ss_ct
eqn_nuprime_ss_ct.nuprime_ss_ct
eqn_y_ss_ct.y_ss_ct
eqn_prof_ss_ct.prof_ss_ct
eqn_phi_ss_ct.phi_ss_ct
eqn_q_ss_ct.q_ss_ct
eqn_theta_ss_ct.theta_ss_ct
eqn_thetaind_ss_ct.thetaind_ss_ct
eqn_thetabar_ss_ct.thetabar_ss_ct
eqn_n_ss_ct.n_ss_ct
eqn1_u_ss_ct
eqn2_u_ss_ct
eqn_ntot_ss_ct.ntot_ss_ct
eqn_ubar_ss_ct.ubar_ss_ct
eqn_tlump_ss_ct.tlump_ss_ct
eqn_mc_ss_ct.p_ss_ct

eqn_pbar.pbar
eqn_pnet.pnet
eqn_hhfoc.c
eqn_hheuler
eqn_hheuler_tc
eqn_hhbc.lambda
eqn_cd.cd
eqn1_u_n.u_n
eqn1_u_u.u_u
eqn_v_n.v_n
eqn_v_n_tc
eqn_v_u.v_u
eqn_v_u_tc
eqn_phibar.phibar
eqn_phibar_tc
eqn_nash.w0,
eqn_wdelta.wdelta,
eqn_wdelta_tc,
eqn_wstar.wstar,
eqn_wstar_tc,
eqn_w.w,
eqn_intra.h
eqn_firmfoc.vbar
eqn_vbarlast
eqn_j_n.j_n
eqn_j_n_tc
eqn_nu.nu
eqn_zbar.zbar
eqn_nuprime.nuprime
eqn_y.y
eqn_prof.prof
eqn_phi.phi
eqn_q.q
eqn_theta.theta
eqn_thetaind.thetaind
eqn_thetabar.thetabar
eqn_n.n
eqn_ntot.ntot
eqn1_u
eqn2_u
eqn_ubar.ubar
eqn_tlump.tlump
eqn_mc.p

eqn_tau_l_ss_ct.tau_l_ss_ct
eqn_tau_p_ss_ct.tau_p_ss_ct
eqn_tau_l.tau_l
eqn_tau_p.tau_p
eqn_diff.delta_tau_l
/;

model edge_b2_taul /
eqn_pbar_ss_ct.pbar_ss_ct
eqn_pnet_ss_ct.pnet_ss_ct
eqn_hhfoc_ss_ct.c_ss_ct
eqn_bc_ss_ct.lambda_ss_ct
eqn_cd_ss_ct.cd_ss_ct
eqn_v_n_ss_ct.v_n_ss_ct
eqn_v_u_ss_ct.v_u_ss_ct
eqn2_u_n_ss_ct.u_n_ss_ct
eqn2_u_u_ss_ct.u_u_ss_ct
eqn_phibar_ss_ct.phibar_ss_ct
eqn_nash_ss_ct.w_ss_ct
eqn_intra_ss_ct.h_ss_ct
eqn_firm_ss_ct.vbar_ss_ct
eqn_j_n_ss_ct.j_n_ss_ct
eqn_nu_ss_ct.nu_ss_ct
eqn_zbar_ss_ct.zbar_ss_ct
eqn_nuprime_ss_ct.nuprime_ss_ct
eqn_y_ss_ct.y_ss_ct
eqn_prof_ss_ct.prof_ss_ct
eqn_phi_ss_ct.phi_ss_ct
eqn_q_ss_ct.q_ss_ct
eqn_theta_ss_ct.theta_ss_ct
eqn_thetaind_ss_ct.thetaind_ss_ct
eqn_thetabar_ss_ct.thetabar_ss_ct
eqn_n_ss_ct.n_ss_ct
eqn1_u_ss_ct
eqn2_u_ss_ct
eqn_ntot_ss_ct.ntot_ss_ct
eqn_ubar_ss_ct.ubar_ss_ct
eqn_tlump_ss_ct.tlump_ss_ct
eqn_mc_ss_ct.p_ss_ct

eqn_pbar.pbar
eqn_pnet.pnet
eqn_hhfoc.c
eqn_hheuler
eqn_hheuler_tc
eqn_hhbc.lambda
eqn_cd.cd
eqn2_u_n.u_n
eqn2_u_u.u_u
eqn_v_n.v_n
eqn_v_n_tc
eqn_v_u.v_u
eqn_v_u_tc
eqn_phibar.phibar
eqn_phibar_tc
eqn_nash.w0,
eqn_wdelta.wdelta,
eqn_wdelta_tc,
eqn_wstar.wstar,
eqn_wstar_tc,
eqn_w.w,
eqn_intra.h
eqn_firmfoc.vbar
eqn_vbarlast
eqn_j_n.j_n
eqn_j_n_tc
eqn_nu.nu
eqn_zbar.zbar
eqn_nuprime.nuprime
eqn_y.y
eqn_prof.prof
eqn_phi.phi
eqn_q.q
eqn_theta.theta
eqn_thetaind.thetaind
eqn_thetabar.thetabar
eqn_n.n
eqn_ntot.ntot
eqn1_u
eqn2_u
eqn_ubar.ubar
eqn_tlump.tlump
eqn_mc.p

eqn_tau_l_ss_ct.tau_l_ss_ct
eqn_tau_p_ss_ct.tau_p_ss_ct
eqn_tau_l.tau_l
eqn_tau_p.tau_p
eqn_diff.delta_tau_l
/;

model edge_b1_taup /
eqn_pbar_ss_ct.pbar_ss_ct
eqn_pnet_ss_ct.pnet_ss_ct
eqn_hhfoc_ss_ct.c_ss_ct
eqn_bc_ss_ct.lambda_ss_ct
eqn_cd_ss_ct.cd_ss_ct
eqn_v_n_ss_ct.v_n_ss_ct
eqn_v_u_ss_ct.v_u_ss_ct
eqn1_u_n_ss_ct.u_n_ss_ct
eqn1_u_u_ss_ct.u_u_ss_ct
eqn_phibar_ss_ct.phibar_ss_ct
eqn_nash_ss_ct.w_ss_ct
eqn_intra_ss_ct.h_ss_ct
eqn_firm_ss_ct.vbar_ss_ct
eqn_j_n_ss_ct.j_n_ss_ct
eqn_nu_ss_ct.nu_ss_ct
eqn_zbar_ss_ct.zbar_ss_ct
eqn_nuprime_ss_ct.nuprime_ss_ct
eqn_y_ss_ct.y_ss_ct
eqn_prof_ss_ct.prof_ss_ct
eqn_phi_ss_ct.phi_ss_ct
eqn_q_ss_ct.q_ss_ct
eqn_theta_ss_ct.theta_ss_ct
eqn_thetaind_ss_ct.thetaind_ss_ct
eqn_thetabar_ss_ct.thetabar_ss_ct
eqn_n_ss_ct.n_ss_ct
eqn1_u_ss_ct
eqn2_u_ss_ct
eqn_ntot_ss_ct.ntot_ss_ct
eqn_ubar_ss_ct.ubar_ss_ct
eqn_tlump_ss_ct.tlump_ss_ct
eqn_mc_ss_ct.p_ss_ct

eqn_pbar.pbar
eqn_pnet.pnet
eqn_hhfoc.c
eqn_hheuler
eqn_hheuler_tc
eqn_hhbc.lambda
eqn_cd.cd
eqn1_u_n.u_n
eqn1_u_u.u_u
eqn_v_n.v_n
eqn_v_n_tc
eqn_v_u.v_u
eqn_v_u_tc
eqn_phibar.phibar
eqn_phibar_tc
eqn_nash.w0,
eqn_wdelta.wdelta,
eqn_wdelta_tc,
eqn_wstar.wstar,
eqn_wstar_tc,
eqn_w.w,
eqn_intra.h
eqn_firmfoc.vbar
eqn_vbarlast
eqn_j_n.j_n
eqn_j_n_tc
eqn_nu.nu
eqn_zbar.zbar
eqn_nuprime.nuprime
eqn_y.y
eqn_prof.prof
eqn_phi.phi
eqn_q.q
eqn_theta.theta
eqn_thetaind.thetaind
eqn_thetabar.thetabar
eqn_n.n
eqn_ntot.ntot
eqn1_u
eqn2_u
eqn_ubar.ubar
eqn_tlump.tlump
eqn_mc.p

eqn_tau_l_ss_ct.tau_l_ss_ct
eqn_tau_p_ss_ct.tau_p_ss_ct
eqn_tau_l.tau_l
eqn_tau_p.tau_p
eqn_diff.delta_tau_p
/;


model edge_b2_taup /
eqn_pbar_ss_ct.pbar_ss_ct
eqn_pnet_ss_ct.pnet_ss_ct
eqn_hhfoc_ss_ct.c_ss_ct
eqn_bc_ss_ct.lambda_ss_ct
eqn_cd_ss_ct.cd_ss_ct
eqn_v_n_ss_ct.v_n_ss_ct
eqn_v_u_ss_ct.v_u_ss_ct
eqn2_u_n_ss_ct.u_n_ss_ct
eqn2_u_u_ss_ct.u_u_ss_ct
eqn_phibar_ss_ct.phibar_ss_ct
eqn_nash_ss_ct.w_ss_ct
eqn_intra_ss_ct.h_ss_ct
eqn_firm_ss_ct.vbar_ss_ct
eqn_j_n_ss_ct.j_n_ss_ct
eqn_nu_ss_ct.nu_ss_ct
eqn_zbar_ss_ct.zbar_ss_ct
eqn_nuprime_ss_ct.nuprime_ss_ct
eqn_y_ss_ct.y_ss_ct
eqn_prof_ss_ct.prof_ss_ct
eqn_phi_ss_ct.phi_ss_ct
eqn_q_ss_ct.q_ss_ct
eqn_theta_ss_ct.theta_ss_ct
eqn_thetaind_ss_ct.thetaind_ss_ct
eqn_thetabar_ss_ct.thetabar_ss_ct
eqn_n_ss_ct.n_ss_ct
eqn1_u_ss_ct
eqn2_u_ss_ct
eqn_ntot_ss_ct.ntot_ss_ct
eqn_ubar_ss_ct.ubar_ss_ct
eqn_tlump_ss_ct.tlump_ss_ct
eqn_mc_ss_ct.p_ss_ct

eqn_pbar.pbar
eqn_pnet.pnet
eqn_hhfoc.c
eqn_hheuler
eqn_hheuler_tc
eqn_hhbc.lambda
eqn_cd.cd
eqn2_u_n.u_n
eqn2_u_u.u_u
eqn_v_n.v_n
eqn_v_n_tc
eqn_v_u.v_u
eqn_v_u_tc
eqn_phibar.phibar
eqn_phibar_tc
eqn_nash.w0,
eqn_wdelta.wdelta,
eqn_wdelta_tc,
eqn_wstar.wstar,
eqn_wstar_tc,
eqn_w.w,
eqn_intra.h
eqn_firmfoc.vbar
eqn_vbarlast
eqn_j_n.j_n
eqn_j_n_tc
eqn_nu.nu
eqn_zbar.zbar
eqn_nuprime.nuprime
eqn_y.y
eqn_prof.prof
eqn_phi.phi
eqn_q.q
eqn_theta.theta
eqn_thetaind.thetaind
eqn_thetabar.thetabar
eqn_n.n
eqn_ntot.ntot
eqn1_u
eqn2_u
eqn_ubar.ubar
eqn_tlump.tlump
eqn_mc.p

eqn_tau_l_ss_ct.tau_l_ss_ct
eqn_tau_p_ss_ct.tau_p_ss_ct
eqn_tau_l.tau_l
eqn_tau_p.tau_p
eqn_diff.delta_tau_p
/;

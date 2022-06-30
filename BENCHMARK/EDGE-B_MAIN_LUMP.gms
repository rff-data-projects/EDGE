* Employment Dynamics in General Equilibrium Benchmark Model - EDGE-BENCHMARK

* Two Sector Model with Search Frictions
* U(c,h) = c^(1-sigma)/(1-sigma) - psi h**(1+chi)/(1+chi);
* Numeraire is price of clean good
* Length of period is one month
* No Capital
* Separation and recruitment costs allowed to vary by sector

* Preexisting taxes
* Sectoral Switching Costs and Sticky Wages

* MAIN_LUMP: Multiple policy simulation, lump-sum rebates


options
        decimals = 5
;

sets
        t "time periods" /t1*t101/
        tp(t) /t1*t101/
        tpi(t) /t1/
        tfirst(t) "First time period"
        tlast(t) "Final time period"
        ind /c,d/
        it "policy iterations" /it1*it2/
;

tfirst(t) = yes$(ord(t) eq 1);
tlast(t) = yes$(ord(t) eq card(t));
alias (ind,ind2)


scalar
        nt "number of periods" /101/
        npi /1/
;


* DEFINE DEEP PARAMETERS
parameters
         beta "time preference" /.996/
         betat(t)
         chi "Frisch Elasticity of Labor Supply" /1/
         eta "Nash bargaining weight" /0.5/
         pi(ind) "Exogenous quit rate"
         gamma "Curvature of matching function" /0.5/
         mu_e(ind) "Emissions factor"
         tau_l0 "Tax on Labor Income" /0.25/
         tau_p0 "Payroll tax" /0.12/
         sigma "Intertemporal Elasticity of Substitution" /1/ 
         sigmac "Elasticity of Substitution b/w goods" /0.75/
         theta1 "Abatement cost scale parameter" /1/
         theta2 "Abatement cost curvature" /2.8/
         xi "Industry switching friction parameter" /1/
         delta(ind,ind2) "Delta kronecker"
         rho "Staggerged wage bargaining parameter"

         c_ss "Total consumption"
         cd_ss(ind) "Consumption by sector"
         alphac(ind) "CES share parameter - consumption"
         gammac "CES scale parameter - consumption"
;


c_ss = 1;
cd_ss('c') = 0.8*c_ss;
cd_ss('d') = 0.2*c_ss;
mu_e('c') = 0.000;
mu_e('d') = .0025;
pi(ind) = 0.035;
betat(t) = beta**(ord(t)-1);
alphac('c') = (0.8)**(1/sigmac) / (0.8**(1/sigmac) + 0.2**(1/sigmac));
alphac('d') = (0.2)**(1/sigmac) / (0.8**(1/sigmac) + 0.2**(1/sigmac));
gammac = 0.8**(1/sigmac) + 0.2**(1/sigmac);

delta(ind2,ind) = 0;
delta('c','c') = 1;
delta('d','d') = 1;

rho = 0;


* CALIBRATE MODEL

$INCLUDE EDGE-B_CAL.gms



* DEFINE EMISSION PRICE PATHS
parameter
p_e0
p_e_ss
p_e(t)
;

p_e0('it1') = 10;
p_e0('it2') = 20;




$INCLUDE EDGE-B_EQN.gms


* REPLICATE BASE CASE WITH P_E = 0;
p_e_ss = 0;
options iterlim = 0;

delta_tau_l.fx = 1;
delta_tau_p.fx = 1;


if(sigma = 1,
solve edge_b1_ss_ct using mcp   ;
);
if(sigma ne 1,
solve edge_b2_ss_ct using mcp   ;
);



* SOLVE STEADY STATE WITH P_E>0
options iterlim = 1000;


p_e_ss = p_e0('it1');

if(sigma = 1,
solve edge_b1_ss_ct using mcp   ;
);
if(sigma ne 1,
solve edge_b2_ss_ct using mcp   ;
);





* STEADY STATE OUTPUT

parameter
l_ss(ind)
l_ss_ct(ind)
perc_l(ind)
perc_ltot
perc_c
perc_n
perc_ntot
perc_nh
perc_nhtot
perc_earntot
e_ss
e_ss_ct
perc_e(it)
ubar_it(it)
;

* TRANSITION OUTPUT

parameters
pve
pve_ref
pv_reduce
ev_ton
tlump_d_ss
crev_ss
net_gross_ss
tlump_d_pv
crev_pv
net_gross_pv
vcost_ss
vcost(t)
pvlump
pvlump_ss
pvlump_diff

perc_pve
pv_reduce_input
pv_reduce_output
emit
e_reduce
e_reduce_out(t)
e_reduce_out_pct(t)

earn
earn_ss_ct
earn_ss
earntot
earntot_ss_ct
earntot_ss
delta_u
delta_n
deltaratio
;

;




* SOLVE FULL MODEL

loop(it,

* Initialize emissions prices
p_e_ss  = p_e0(it);
p_e(t) = p_e_ss;

* Initialize first transition with steady state output

if(ord(it) = 1,
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
u.l(ind,t) = u_ss_ct.l(ind);
ntot.l(t) = ntot_ss_ct.l;
ubar.l(t) = ubar_ss_ct.l;
tlump.l(t) = tlump_ss_ct.l;
;
);





if(sigma = 1,
solve edge_b1_lump using mcp;
);
if(sigma ne 1,
solve edge_b2_lump using mcp;
);


* CALCULATE STEADY STATE OUTPUT

l_ss(ind) = h_ss.l(ind)*n_ss.l(ind)*(1-vbar_ss.l(ind));

l_ss_ct(ind) = h_ss_ct.l(ind)*n_ss_ct.l(ind)*(1-vbar_ss_ct.l(ind));
perc_l(ind) = 100*(l_ss_ct(ind)-l_ss(ind))/l_ss(ind);
perc_ltot = 100*(sum(ind,l_ss_ct(ind)) - sum(ind,l_ss(ind)))/sum(ind,l_ss(ind));
perc_c = 100*(c_ss_ct.l-c_ss)/c_ss;
perc_n(ind) = 100*(n_ss_ct.l(ind)-n_ss.l(ind))/n_ss.l(ind);
perc_ntot = 100*(sum(ind,n_ss_ct.l(ind)) - sum(ind,n_ss.l(ind)))/sum(ind,n_ss.l(ind));
perc_nh(ind) = 100*(n_ss_ct.l(ind)*h_ss_ct.l(ind)-n_ss.l(ind)*h_ss.l(ind))/(n_ss.l(ind)*h_ss.l(ind));
perc_nhtot = 100*(sum(ind,n_ss_ct.l(ind)*h_ss_ct.l(ind)) - sum(ind,n_ss.l(ind)*h_ss.l(ind)))/sum(ind,n_ss.l(ind)*h_ss.l(ind));
perc_earntot = 100*(sum(ind,(1-tau_l_ss_ct.l)*w_ss_ct.l(ind)*n_ss_ct.l(ind)*h_ss_ct.l(ind))/pbar_ss_ct.l - sum(ind,(1-tau_l0)*w_ss.l(ind)*n_ss.l(ind)*h_ss.l(ind))/pbar_ss)/(sum(ind,(1-tau_l0)*w_ss.l(ind)*n_ss.l(ind)*h_ss.l(ind))/pbar_ss);

e_ss = sum(ind,mu_e(ind)*y_ss.l(ind));
e_ss_ct = sum(ind,(1-nu_ss_ct.l(ind))*mu_e(ind)*y_ss_ct.l(ind));
perc_e(it) = 100*(e_ss_ct-e_ss)/e_ss;

ubar_it(it) = ubar_ss_ct.l;


display
perc_l
perc_ltot
perc_c
perc_n
perc_ntot
perc_nh
perc_nhtot
perc_earntot
ubar_ss_ct.l
perc_e
zbar_ss_ct.l
;


* CALCULATE TRANSITION OUTPUT




pve = sum(t,sum(ind,((1-nu.l(ind,t))*mu_e(ind)*y.l(ind,t)))*betat(t)) + (sum(ind,((1-nu_ss_ct.l(ind))*mu_e(ind)*y_ss_ct.l(ind)))*beta**nt) / (1-beta);
pve_ref = sum(ind,(mu_e(ind)*y_ss.l(ind))) / (1-beta);
pv_reduce = pve_ref - pve;


tlump_d_ss = tlump_ss_ct.l - tlump_ss;
crev_ss = sum(ind,(1-nu_ss_ct.l(ind))*p_e_ss*mu_e(ind)*y_ss_ct.l(ind)/pbar_ss_ct.l);
net_gross_ss = tlump_d_ss / crev_ss;

tlump_d_pv  = sum(t,tlump.l(t)*betat(t)) + (tlump_ss_ct.l*beta**nt)/(1-beta) - tlump_ss/(1-beta);
crev_pv = sum(t,sum(ind,((1-nu.l(ind,t))*p_e(t)*mu_e(ind)*y.l(ind,t))/pbar.l(t) )*betat(t)) + (sum(ind,((1-nu_ss_ct.l(ind))*p_e_ss*mu_e(ind)*y_ss_ct.l(ind))/pbar_ss_ct.l)*beta**nt)/(1-beta);
net_gross_pv = tlump_d_pv / crev_pv;


perc_pve = 100*pv_reduce/pve_ref;
pv_reduce_input =  pve_ref - sum(t,sum(ind,((1-nu.l(ind,t))*mu_e(ind)*y_ss.l(ind)))*betat(t)) - (sum(ind,((1-nu_ss_ct.l(ind))*mu_e(ind)*y_ss.l(ind)))*beta**nt) / (1-beta);
pv_reduce_output =  sum(t,sum(ind,((1-nu.l(ind,t))*mu_e(ind)*(y_ss.l(ind) - y.l(ind,t))))*betat(t)) + (sum(ind,((1-nu_ss_ct.l(ind))*mu_e(ind)*(y_ss.l(ind)-y_ss_ct.l(ind))))*beta**nt) / (1-beta);

emit(t) = (1-nu.l('d',t))*mu_e('d')*y.l('d',t);

e_reduce(t) = e_ss - emit(t);
e_reduce_out(t) = e_ss - mu_e('d')*y.l('d',t);
e_reduce_out_pct(t) = e_reduce_out(t)/e_reduce(t);




display
pv_reduce
pv_reduce_input
pv_reduce_output
net_gross_ss
net_gross_pv
tlump_d_ss
tlump_d_pv
vbar_ss.l
vbar_ss_ct.l
perc_pve
ubar.l
;




earn_ss(ind) = (1-tau_l0)*w_ss.l(ind)*h_ss.l(ind)/pbar_ss;
earn_ss_ct(ind) = (1-tau_l_ss_ct.l)*w_ss_ct.l(ind)*h_ss_ct.l(ind)/pbar_ss_ct.l;
earn(ind,t) = (1-tau_l.l(t))*w.l(ind,t)*h.l(ind,t)/pbar.l(t);
earntot_ss(ind) = (1-tau_l0)*w_ss.l(ind)*h_ss.l(ind)*n_ss.l(ind)/pbar_ss;
earntot_ss_ct(ind) = (1-tau_l_ss_ct.l)*w_ss_ct.l(ind)*h_ss_ct.l(ind)*n_ss_ct.l(ind)/pbar_ss_ct.l;
earntot(ind,t) = (1-tau_l.l(t))*w.l(ind,t)*h.l(ind,t)*n.l(ind,t)/pbar.l(t);


delta_u(it) = ubar_ss_ct.l - ubar_ss.l;
delta_n(it) = (n_ss.l('d') - n_ss_ct.l('d')) + (n_ss_ct.l('c')-n_ss.l('c'));
deltaratio(it) = delta_n(it) / delta_u(it);


display
delta_u
delta_n
deltaratio
ubar_it
;

);

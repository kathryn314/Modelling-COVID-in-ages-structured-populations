function sir_simulation

beta = zeros(3); % rate of infection (or transmission parameter) which is same for both strains

beta(1,1) = 0.8;
beta(1,2) = 0.02;
beta(1,3) = 0.01;
beta(2,1) = 0.02;
beta(2,2) = 0.8;
beta(2,3) = 0.03;
beta(3,1) = 0.01;
beta(3,2) = 0.03;
beta(3,3) = 0.8;

gamma = 0.0476; % rate of recovery
mu = 0.000055; % rate of host natural death

nu = zeros(4,1); % aging rate

nu(2) = 0.00001; % 1/aging rate = average time somebody spends in the source compartment (e.g. ages 1-20 then ageing rate is 1/20 years)
nu(3) = 0.00001;

N = [100 ,100 ,100]; % Total population N = S + I + R

p = [0.9, 0.5, 0.9]; % vaccine coverage

xi = 0.00014; % rate of loss of natural immunity (assumption that this loss stays the same at all age groups)

xi_v = 0.00014; % rate of loss of vaccine immunity

theta = 0.6; % primary reduced susceptibility

tau = 0.2; % vaccinated reduced susceptibility

v = 0; % primary reduced infection period

n = 0; % vaccinated reduced infection period

B = reshape(beta,[1,9]);

pars = [B,gamma, mu, nu(2), nu(3), N, p, xi, xi_v, theta, tau, v, n];

V_0 = zeros(1,3); % initial number of vaccinated and still immune to strain 1

S_0 = [90 100 100]; % initial number of susceptibles in population

I_0 = [10 0 0 ; 0 0 0]; % initial number of infectives in population

J_0 = zeros(2,3); % initial number of secondary infected in population

R_0 = zeros(3); % initial number of recovered individauls in population

IV2_0 = zeros(1,3); % initial number vaccinated and infected by strain 2

dt = 1; % time interval (days)

t_initial = 0;
t_final = 10000;
tspan = [t_initial : dt: t_final];

y_0 = [V_0(1,:),S_0(1,:), I_0(1,:),I_0(2,:),J_0(1,:), J_0(2,:), R_0(1,:),R_0(2,:),R_0(3,:), IV2_0(1,:)];

[t,y] = ode45(@sir_model, tspan, y_0,[], pars);

hold on
plot(t,y(:,10)); % plot of prevalence of infection in age group 1
plot(t,y(:,11)); % plot of prevalence of infection in age group 2
plot(t,y(:,12)); % plot of prevalence of infection in age group 3

legend({'Infection Age 1','Infection Age 2', 'Infection Age 3'},'Location','northeast')

% plot(t,y(:,1)); % plot of prevalence of vaccination in age group 1
% plot(t,y(:,2)); % plot of prevalence of vaccination in age group 2
% plot(t,y(:,3)); % plot of prevalence of vaccination in age group 3

%legend({'Vaccination Age 1','Vaccination Age 2', 'Vaccination Age 3'},'Location','northeast')

% plot(t,y(:,10)); % plot of prevalence of infection in age group 1
% plot(t,y(:,11)); % plot of prevalence of infection in age group 2
% plot(t,y(:,12)); % plot of prevalence of infection in age group 3

% legend({'Infection Age 1','Infection Age 2', 'Infection Age 3'},'Location','northeast')
end

function f = sir_model(t,y,pars)

f=zeros(15,1);

beta_11 = pars(1);
beta_12 = pars(2);
beta_13 = pars(3);
beta_21 = pars(4);
beta_22 = pars(5);
beta_23 = pars(6);
beta_31 = pars(7);
beta_32 = pars(8);
beta_33 = pars(9);

gamma = pars(10);
mu = pars(11);
nu_1 = pars(12);
nu_2 = pars(13);
N_1 = pars(14);
N_2 = pars(15);
N_3 = pars(16);
p_1 = pars(17);
p_2 = pars(18);
p_3 = pars(19);
xi = pars(20);
xi_v = pars(21);
theta = pars(22);
tau = pars(23);
v = pars(24);
n = pars(25);

S_1 = y(1);
S_2 = y(2);
S_3 = y(3);
I1_1 = y(4);
I1_2 = y(5);
I1_3 = y(6);
I2_1 = y(7);
I2_2 = y(8);
I2_3 = y(9);
J1_1 = y(10);
J1_2 = y(11);
J1_3 = y(12);
J2_1 = y(13);
J2_2 = y(14);
J2_3 = y(15);
R1_1 = y(16);
R1_2 = y(17);
R1_3 = y(18);
R2_1 = y(19);
R2_2 = y(20);
R2_3 = y(21);
R_1 = y(22);
R_2 = y(23);
R_3 = y(24);
V_1 = y(25);
V_2 = y(26);
V_3 = y(27);
IV2_1 = y(28);
IV2_2 = y(29);
IV2_3 = y(30);

lambda1_1 = beta_11*(I1_1+J1_1)/N_1 + beta_12*(I1_2+J1_2)/N_2 + beta_13*(I1_3+J1_3)/N_3;
lambda1_2 = beta_21*(I1_1+J1_1)/N_1 + beta_22*(I1_2+J1_2)/N_2 + beta_23*(I1_3+J1_3)/N_3;
lambda1_3 = beta_31*(I1_1+J1_1)/N_1 + beta_32*(I1_2+J1_2)/N_2 + beta_33*(I1_3+J1_3)/N_3;

lambda2_1 = beta_11*(I2_1+J2_1 +IV2_1)/N_1 + beta_12*(I2_2+J2_2+IV2_2)/N_2 + beta_13*(I2_3+J2_3+IV2_3)/N_3;
lambda2_2 = beta_21*(I2_1+J2_1+IV2_1)/N_1 + beta_22*(I2_2+J2_2+IV2_2)/N_2 + beta_23*(I2_3+J2_3+IV2_3)/N_3;
lambda2_3 = beta_31*(I2_1+J2_1+IV2_1)/N_1 + beta_32*(I2_2+J2_2+IV2_2)/N_2 + beta_33*(I2_3+J2_3+IV2_3)/N_3;

f(1) = (1-p_1)*mu - (lambda1_1+lambda2_1+mu)*S_1 - nu_1*S_1 + xi_v*V_1 + xi*(R1_1+R2_1+R_1);
f(2) = (1-p_2)*mu - (lambda1_2+lambda2_2+mu)*S_2 + nu_1*S_1 - nu_2*S_2 + xi_v*V_2 + xi*(R1_2+R2_2+R_2);
f(3) = (1-p_3)*mu - (lambda1_3+lambda2_3+mu)*S_3 + nu_2*S_2 + xi_v*V_3 + xi*(R1_3+R2_3+R_3);

f(4) = lambda1_1*S_1 - (gamma + mu + nu_1)*I1_1;
f(5) = lambda1_2*S_2 - (gamma + mu +nu_2)*I1_2 + nu_1*I1_1;
f(6) = lambda1_3*S_3 - (gamma+mu)*I1_3 + nu_2*I1_2;

f(7) = lambda2_1*S_1 - (gamma + mu + nu_1)*I2_1;
f(8) = lambda2_2*S_2 - (gamma + mu +nu_2)*I2_2 + nu_1*I2_1;
f(9) = lambda2_3*S_3 - (gamma+mu)*I2_3 + nu_2*I2_2;

f(10) = lambda1_1*(1-theta)*R2_1 - (gamma/(1-v) + mu + nu_1)*J1_1;
f(11) = lambda1_2*(1-theta)*R2_2 - (gamma/(1-v) + mu +nu_2)*J1_2 + nu_1*J1_1;
f(12) = lambda1_3*(1-theta)*R2_3 - (gamma/(1-v)+ mu)*J1_3 + nu_2*J1_2;

f(13) = lambda2_1*(1-theta)*R1_1 - (gamma/(1-v) + mu + nu_1)*J2_1;
f(14) = lambda2_2*(1-theta)*R1_2 - (gamma/(1-v) + mu +nu_2)*J2_2 + nu_1*J2_1;
f(15) = lambda2_3*(1-theta)*R1_3 - (gamma/(1-v)+mu)*J2_3 + nu_2*J2_2;

f(16) = gamma*I1_1 - lambda2_1*(1-theta)*R1_1 - (mu + xi)*R1_1 - nu_1*R1_1;
f(17) = gamma*I1_2 - lambda2_2*(1-theta)*R1_2 - (mu + xi)*R1_2 - nu_2*R1_2 +nu_1*R1_1;
f(18) = gamma*I1_3 - lambda2_3*(1-theta)*R1_3 - (mu + xi)*R1_3 + nu_2*R1_2;

f(19) = gamma*I2_1 - lambda1_1*(1-theta)*R2_1 - (mu + xi)*R2_1 - nu_1*R2_1;
f(20) = gamma*I2_2 - lambda1_2*(1-theta)*R2_2 - (mu + xi)*R2_2 - nu_2*R2_2 +nu_1*R2_1;
f(21) = gamma*I2_3 - lambda1_3*(1-theta)*R2_3 - (mu + xi)*R2_3 + nu_2*R2_2;

f(22) = gamma*((J1_1+J2_1)/(1-v) + IV2_1/(1-n)) - (mu+nu_1+xi)*R_1;
f(23) = gamma*((J1_2+J2_2)/(1-v) + IV2_2/(1-n)) - (mu+nu_2+xi)*R_2 + nu_1*R_1;
f(24) = gamma*((J1_3+J2_3)/(1-v) + IV2_3/(1-n)) - (mu+xi)*R_3 + nu_2*R_2;

f(25) = (p_1 - V_1)*mu - nu_1*V_1 - xi_v*V_1 - lambda2_1*(1-tau)*V_1;
f(26) = (p_2 - V_2)*mu + nu_1*V_1 - nu_2*V_2 - xi_v*V_2- lambda2_2*(1-tau)*V_2;
f(27) = (p_3 - V_3)*mu + nu_2*V_2 - xi_v*V_3 - lambda2_3*(1-tau)*V_3;

f(28) = lambda2_1*(1-tau)*V_1 - (gamma/(1-n)+mu)*IV2_1 -nu_1*IV2_1;
f(29) = lambda2_2*(1-tau)*V_2 - (gamma/(1-n)+mu)*IV2_2 -nu_2*IV2_2 + nu_1*IV2_1;
f(30) = lambda2_3*(1-tau)*V_3 - (gamma/(1-n)+mu)*IV2_3 +nu_2*IV2_2;

end
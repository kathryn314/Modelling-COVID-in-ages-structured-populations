function sir_simulation

beta = zeros(3); % rate of infection (or transmission parameter)

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
sigma = 0.2; % incubation rate

nu = zeros(4,1); % aging rate

nu(2) = 0.00001;
nu(3) = 0.00001;

N = [100 ,10 ,10]; % Total population N = S + I + R

p = [0.9, 0.5, 0.9]; % vaccine coverage

xi = 0.00014; % rate of loss of natural immunity (assumption that this loss stays the same at all age groups)

xi_v = 0.00014; % rate of loss of vaccine immunity

B = reshape(beta,[1,9]);

pars = [B,gamma, mu, sigma, nu(2), nu(3), N, p, xi, xi_v];

V_0 = zeros(3,1); % initial number of vaccinated in population

S_0 = [90 10 10]; % initial number of susceptibles in population

E_0 = zeros(3,1); % initial number latent in population

I_0 = [10 0 0]; % initial number of infectives in population

R_0 = zeros(3,1); % initial number of recovered individauls in population

dt = 1; % time interval (1/4 of a day)

t_initial = 0;
t_final = 1000;
tspan = [t_initial : dt: t_final];

y_0 = [V_0(1),V_0(2),V_0(3), S_0(1),S_0(2),S_0(3), E_0(1),E_0(2),E_0(3), I_0(1),I_0(2),I_0(3), R_0(1),R_0(2),R_0(3)];

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
sigma = pars(12);
nu_1 = pars(13);
nu_2 = pars(14);
N_1 = pars(15);
N_2 = pars(16);
N_3 = pars(17);
p_1 = pars(18);
p_2 = pars(19);
p_3 = pars(20);
xi = pars(21);
xi_v = pars(22);

V_1 = y(1);
V_2 = y(2);
V_3 = y(3);
S_1 = y(4);
S_2 = y(5);
S_3 = y(6);
E_1 = y(7);
E_2 = y(8);
E_3 = y(9);
I_1 = y(10);
I_2 = y(11);
I_3 = y(12);
R_1 = y(13);
R_2 = y(14);
R_3 = y(15);

lambda_1 = beta_11*I_1/N_1 + beta_12*I_2/N_2 + beta_13*I_3/N_3;
lambda_2 = beta_21*I_1/N_1 + beta_22*I_2/N_2 + beta_23*I_3/N_3;
lambda_3 = beta_31*I_1/N_1 + beta_32*I_2/N_2 + beta_33*I_3/N_3;

f(1) = (p_1 - V_1)*mu - nu_1*V_1 - xi_v*V_1;
f(2) = (p_2 - V_2)*mu + nu_1*V_1 - nu_2*V_2 - xi_v*V_2;
f(3) = (p_3 - V_3)*mu + nu_2*V_2 - xi_v*V_2;

f(4) = (1-p_1)*mu - (lambda_1+mu)*S_1 - nu_1*S_1 + xi_v*V_1;
f(5) = (1-p_2)*mu- (lambda_2+mu)*S_2 + nu_1*S_1 - nu_2*S_2 + xi_v*V_2;
f(6) = (1-p_3)*mu- (lambda_3+mu)*S_3 + nu_2*S_2 + xi_v*V_3;

f(7) = lambda_1*S_1 - (sigma + mu + nu_1)*E_1;
f(8) = lambda_2*S_2 - (sigma + mu +nu_2)*E_2 + nu_1*E_1;
f(9) = lambda_3*S_3 - (sigma+mu)*E_3 + nu_2*E_2;

f(10) = sigma*E_1 - (gamma+mu+nu_1)*I_1;
f(11) = sigma*E_2 - (gamma+mu+nu_2)*I_2+nu_1*I_1;
f(12) = sigma*E_3 - (gamma+mu)*I_3+nu_2*I_2;

f(13) = gamma*I_1 - (mu+nu_1+xi)*R_1;
f(14) = gamma*I_2 - (mu+nu_2+xi)*R_2+nu_1*R_1;
f(15) = gamma*I_3 - (mu+xi)*R_3+nu_2*R_2;
end
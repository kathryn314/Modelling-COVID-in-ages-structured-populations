function sir_simulation

mu = 0.000055; % natural birth and death rate
beta = 0.81; % rate of transmission
gamma = 0.0476; % rate of recovery
N = 1; % total population N = S + I + R - keep N = 1 since we have adapted our model for this to be the case
p = 1; % vaccine coverage
sigma = 0.00014; % rate of loss of natural immunity 
sigma_v = 0.00014; % rate of loss of vaccine immunity
theta = 0.6; % primary reduced susceptibility
v = 0; % primary reduced infection period
tau = 0.2; % vaccinated reduced susceptibility
nu = 0; % vaccinated reduced infection period

pars = [mu, beta, gamma, sigma, sigma_v, N, p, theta, v, tau, nu];

V_0 = p*mu/(mu + sigma_v); % initial number vaccinated and still immune to strain 1
S_0 = (gamma+mu)/beta % initial number of susceptibles in population
I_0 = [(mu+sigma)*(1-S_0-V_0)/(gamma+mu+sigma) 10^(-6)] % initial number of infectives in population
J_0 = [0 0]; % initial number of secondary infected in population
R_0 = [gamma*(1-S_0-V_0)/(gamma+mu+sigma) 0 0] % initial number of recovered individauls (and therefore immune) in population
IV2_0 = 0; % initial number vaccinated and infected by strain 2

dt = 1; % time interval (days)

t_initial = 0;
t_final = 1500;
tspan = [t_initial : dt: t_final];

y_0 = [S_0;I_0(1);I_0(2); J_0(1); J_0(2) ;R_0(1); R_0(2); R_0(3); V_0; IV2_0];

[t,y] = ode45(@sir_model, tspan, y_0,[], pars);

hold on
plot(t,log(y(:,2))); % plot of prevalence of first infection by strain 1 over time
% plot(t,y(:,3)); % plot of prevalence of first infection by strain 1 over time
% plot(t,y(:,4)); % plot of prevalence of second infection by strain 2 over time
% plot(t,y(:,5)); % plot of prevalence of second infection by strain 1 over time
% plot(t,y(:,1)); % plot of prevalence of susceptability over time
% plot(t,y(:,6)); % plot of prevalence of recovery from strain 1 over time
% plot(t,y(:,7)); % plot of prevalence of recovery from strain 2 over time
% plot(t,y(:,8)); % plot of prevalence of recovery from both strains over time

legend({'Strain 1','Strain 2', 'Strain 12', 'Strain 21'},'Location','northeast')


end

function f = sir_model(t,y,pars)

f=zeros(10,1);
mu = pars(1);
beta = pars(2);
gamma = pars(3);
sigma = pars(4);
sigma_v = pars(5);
N = pars(6);
p = pars(7);
theta = pars(8); 
v = pars(9);
tau = pars(10);
nu = pars(11);


S = y(1);
I_1 = y(2);
I_2 = y(3);
J_1 = y(4);
J_2 = y(5);
R_1 = y(6);
R_2 = y(7);
R = y(8);
V = y(9);
IV2 = y(10);

f(1) = mu*(1-p-S) - (I_1+J_1+I_2+J_2+IV2)*beta*S+sigma*(R_1+R_2+R) +sigma_v*V;
f(2) = beta*S*(I_1+J_1)-(gamma+mu)*I_1;
f(3) = beta*S*(I_2+J_2+IV2) - (gamma+mu)*I_2;
f(4) = beta*(1-theta)*R_2*(I_1+J_1) - (gamma/(1-v)+mu)*J_1;
f(5) = beta*(1-theta)*R_1*(I_2+J_2+IV2) - (gamma/(1-v)+mu)*J_2;
f(6) = gamma*I_1 - beta*(1-theta)*R_1*(I_2+J_2+IV2) - (mu+sigma)*R_1;
f(7) = gamma*I_2 - beta*(1-theta)*R_2*(I_1+J_1) - (mu+sigma)*R_2;
f(8) = gamma*((J_1+J_2)/(1-v)+IV2/(1-nu))-(mu+sigma)*R;
f(9) = mu*(p-V) - beta*(1-tau)*V*(I_2+J_2+IV2) - sigma_v*V;
f(10) = beta*(1-tau)*V*(I_2+J_2+IV2) - (gamma/(1-nu)+mu)*IV2;

end

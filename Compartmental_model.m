function sir_simulation

beta = 1; % rate of transmission between susceptible and infected individual
gamma = 0.05; % rate of recovery
delta = 0; % rate of loss of immunity
mu = 0.05; % rate of host natural death
alpha = 0; % virulence
N = 1; % Total population N = S + I + R

pars = [alpha, beta, gamma, delta, mu, N];

I_0 = 0.001; % initial number of infectives in population
S_0 = 0.999; % initial number of susceptibles in population
R_0 = 0; %initial number of recovered individauls in population

dt = 0.25; % time interval (1/4 of a day)

fprintf('Value of parameter R_0 is %.2f' , beta/(gamma + alpha + mu))

t_initial = 0;
t_final = 100;
tspan = [t_initial : dt: t_final];

y_0 = [S_0;I_0;R_0];

[t,y] = ode45(@sir_model, tspan, y_0,[], pars);

hold on
plot(t,y(:,2)); % plot of prevalence of infection over time
plot(t,y(:,1)); % plot of prevalence of susceptability over time
plot(t,y(:,3)); % plot of prevalence of recovery over time
end

function f = sir_model(t,y,pars)

f=zeros(3,1);
alpha = pars(1);
beta = pars(2);
gamma = pars(3);
delta = pars(4);
mu = pars(5);
N = pars(6);

S = y(1);
I = y(2);
R = y(3);

f(1) = - mu*S - beta*S*I/N + delta*R;
f(2) = beta*S*I/N - (gamma + alpha + mu)*I;
f(3) = gamma*I - (delta + mu)*R;
end

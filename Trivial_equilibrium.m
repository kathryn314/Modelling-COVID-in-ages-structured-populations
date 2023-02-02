function sir_simulation

beta = [29*10^(-4) 29*10^(-4); 29*10^(-4) 29*10^(-4)]; % rate of infection (or transmission parameter)
gamma = [26 100]; % rate of recovery
delta = 0.02; % mortality rate of uninfected hosts
alpha = [0 0]; % virulence
N = 10^5; % total population N = S + I + R
p = 0.5; % proportion vaccinated
b = delta*N; % flux arrival of new susceptible hosts

pars = [alpha(1), alpha(2), beta(1,1), beta(1,2), beta(2,1), beta(2,2), gamma(1), gamma(2), delta, b, N, p];

RN_0 = (beta(1,1)/(delta+alpha(1)+gamma(1)))*b/delta; 
RV_0 = (beta(2,2)/(delta+alpha(2)+gamma(2)))*b/delta;
R_0 = (1-p)*RN_0+p*RV_0; 

pc = 1- (1-RV_0)/(RN_0-RV_0);

I_0 = [(delta/beta(1,1))*(RN_0 - 1) 0]; % initial proportion of infectives in population
S_0 = [b/(delta*RN_0) 0]; % initial proportion of susceptibles in population
Re_0 = (gamma(1)/beta(1,1))*(RN_0 - 1); % initial proportion of recovered individauls in population

dt = 0.25; % time interval

fprintf('Value of parameter R_0 is %.2f' , R_0)

t_initial = 0;
t_final = 100;
tspan = [t_initial : dt: t_final];

y_0 = [S_0(1);S_0(2);I_0(1); I_0(2);Re_0];

[t,y] = ode45(@sir_model, tspan, y_0,[], pars);

hold on
plot(t,log10(y(:,3)+y(:,4))); % plot of prevalence of infection over time
% plot(t,y(:,1)+y(:,2)); % plot of prevalence of susceptability over time
% plot(t,y(:,5)); % plot of prevalence of recovery over time

% plot(t,y(:,3))
% plot(t,y(:,4))

end

function f = sir_model(t,y,pars)

f=zeros(5,1);

alpha_N = pars(1);
alpha_V = pars(2);
beta_NN = pars(3);
beta_NV = pars(4);
beta_VN = pars(5);
beta_VV = pars(6);
gamma_N = pars(7);
gamma_V = pars(8);
delta = pars(9);
lambda = pars(10);
N = pars(11);
p = pars(12);

S_N = y(1);
S_V = y(2);
I_N = y(3);
I_V = y(4);
Re = y(5);

L_N = beta_NN*I_N + beta_VN*I_V;
L_V = beta_VV*I_V + beta_NV*I_N;

f(1) = lambda*(1-p) - (delta + L_N)*S_N;
f(2) = lambda*p - (delta + L_V)*S_V;
f(3) = L_N*S_N - (delta + alpha_N + gamma_N)*I_N;
f(4) = L_V*S_V - (delta + alpha_V + gamma_V)*I_V;
f(5) = gamma_N*I_N+gamma_V*I_V- delta*Re;
end

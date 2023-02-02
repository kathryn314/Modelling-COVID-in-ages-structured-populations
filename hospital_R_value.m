gamma = 0.0476
    tau = 0.5
    beta = 0.81

for nu = [0.01 0.0238 0.1 1]
    
x = linspace(0,1)

y = (beta./(gamma+x)).*(1+(tau*x)/nu)

hold on
plot(x,y)

legend({'\eta = 0.01','\eta = 0.0238', '\eta = 0.1', '\eta = 1'},'Location','northeast')

xlabel('\xi')
ylabel('Basic reproductive ratio')

print -depsc hospvarying.eps


end 


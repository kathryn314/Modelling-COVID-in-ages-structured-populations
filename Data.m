y = [ 1 2 5 4 8 12 5 22 39 56 56 50 81 67 57 147 266 405 481 479 363 442 609 768 999 1050 1263 1191 1374 2320 2371 2692]'

hold on

plot(1:length(y),y)

x=0:0.01:32
Y = 13.7*exp(0.165*x)

plot(x,Y)

xlabel('Time ( days )')
ylabel('Daily cases')

legend({'Data','y = 13.7 exp(0.165x)'},'Location','northeast')


print -depsc Growthrate.eps
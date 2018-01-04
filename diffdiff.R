#Dests = array(5)
#for (j in 1:5){

set.seed(3)
N = 12500;
dt = .001;
D=1;

xi1 = rnorm(N);
xi2= rnorm(N);
xi3= rnorm(N);


xx = cumsum(xi1*sqrt(2*D*dt))
yy = cumsum(xi2*sqrt(2*D*dt))
tvals = 0:(N-1)*dt;



tau=2;
sigma=1;

Yvals = array(N);
Yvals[1]=0;

for (j in 2:N){
  Yvals[j]= Yvals[j-1] - dt*Yvals[j-1]/tau + sigma*sqrt(dt)*xi3[j];
}

Dvals = Yvals^2;

xxx = cumsum(xi1*sqrt(2*Dvals*dt))
yyy = cumsum(xi2*sqrt(2*Dvals*dt))

plot(x, y,'l',xlim=c(min(xx,x),max(xx,x)), ylim=c(min(yy,y),max(yy,y)))
lines(xx,yy,col=5)



# #

 M = 500;
 msdvals1 = array(M);
 msdvals2 = array(M);


 for (i in 1:M) {
   dix = diff(xx, lag = i)
   diy = diff(yy,lag=i)
   disp1 = sqrt(dix^2+diy^2)
   msd1 = mean(disp1^2)
   
   dixx = diff(xxx, lag = i)
   diyy = diff(yyy,lag=i)
   disp2 = sqrt(dixx^2+diyy^2)
   msd2 = mean(disp2^2)

   msdvals1[i] = msd1
   msdvals2[i] = msd2
 }

 lags = dt * (1:M)
 
 plot(1:M, msdvals1, type='l',lty=2,lwd=5)
 lines(1:M, msdvals2, type='l',lty=3,lwd=5,col=5)
 
 linfit = lm(msdvals1 ~ lags)
 Dest = linfit$coefficients[2] / 4
 print(Dest)
 
 dev.off()
 disp11 = hist(abs(diff(xx, lag = 1) + diff(yy, lag = 1)))
 disp12 = hist(abs(diff(xx, lag = 2) + diff(yy, lag = 2)))
 disp13 = hist(abs(diff(xx, lag = 5) + diff(yy, lag = 5)))
 disp14 = hist(abs(diff(xx, lag = 10) + diff(yy, lag = 10)))
 

 disp21 = hist(abs(diff(xxx, lag = 1) + diff(yyy, lag = 1)))
 disp22 = hist(abs(diff(xxx, lag = 2) + diff(yyy, lag = 2)))
 disp23 = hist(abs(diff(xxx, lag = 5) + diff(yyy, lag = 5)))
 disp24 = hist(abs(diff(xxx, lag = 10) + diff(yyy, lag = 10)))
 
 d11m = c(-rev(disp11$mids),disp11$mids)
 d11d = c(rev(disp11$density),disp11$density)/2
 
 
 plot(disp11$mids,disp11$density,type='l', lwd=6,log='y',xlim=c(0,.4),ylim=c(1e-1,2e1))
 lines(disp12$mids,disp12$density,lwd=3)
 lines(disp13$mids,disp13$density,lwd=1)
 lines(disp14$mids,disp14$density,lwd=0.5)
 
  
 lines(disp21$mids,disp21$density,lwd=6,col=5)
 lines(disp22$mids,disp22$density,lwd=3,col=5)
 lines(disp23$mids,disp23$density,lwd=1,col=5)
 lines(disp24$mids,disp24$density,lwd=0.5,col=5)
 
# #Dests[j]=Dest

# #hist(disp)

# #}

 #h1 = hist(abs(diff(x, lag = 1)))
 #h2 = hist(abs(diff(x, lag = 10)))
 #h3 = hist(abs(diff(x, lag = 25)))
 #h4 = hist(abs(diff(x, lag = 100)))
 #h5 = hist(abs(diff(x, lag = 250)))

# plot(
#   h1$mids,
#   h1$density,
#   type = 'l',
#   lty = 1,
#   lwd = 2.5,
#   xlim = c(0, 2),
#   log = 'y'
# )
# lines(h2$mids, h2$density, lwd = 2)
# lines(h3$mids, h3$density, lwd = 1.5)
# lines(h4$mids, h4$density, lwd = 1)
# lines(h5$mids, h5$density, lwd = 0.5)


#print('---')
#print(mean(Dests))
dt = .001; # step size
T = 1; # max time
D = 1; # diff coefficient
b = .6;#.01; # decay rate
L = 1; # length of domain
a=5;

Nt = round(T/dt) # number of time stpes
dx = .05; 
Nx = round(L/dx) # number of spatial points
x = seq(0, L, length.out=Nx+1)       # Mesh points in space
u = numeric(Nx+1) # empty vector of solutions



# exact solution to PDE, easily obtained via Fourier 
exact <- function(x, t){
  k = pi/L;
  return (exp(-(D*k^2 + b)*t)*sin(k*x))
}

u = exact(x,0) # set initial condition to sin(kx)

## reaction just decay f(u,x,t) = -bu
f <- function(u, x,t){
 return (a+-b*u)
}

# forward euler, u_{i+1} = u_i + dt*f(u_i)
rxn_step = function(u,dt){ 
  return (u + dt*f(u))
}

rxn_step_AB = function(u,uprev,dt){ 
    return (u + 0.5*dt*(3*f(u)-f(uprev)))
}


# Construct tridiagonal matrix
tridiag <- function(upper, lower, main){
    out <- matrix(0,length(main),length(main))
    diag(out) <- main
    indx <- seq.int(length(upper))
    out[cbind(indx+1,indx)] <- lower
    out[cbind(indx,indx+1)] <- upper
    return(out)
} 


# Crank Nicolson, b=0 is just diffusion, b!=0 includes rxns
diffusion_step = function(u,bval,dt){
    mu = 0.5*D*dt/(dx^2); # CFL number
    diagonal = rep(1+2*mu, Nx+1)
    lower = rep(-mu,Nx)
    upper = rep(-mu,Nx)
    
    # BCs
    diagonal[1] = 1;
    upper[1] = 0;
    diagonal[Nx+1] = 1;
    lower[Nx] = 0;
    A<- tridiag(upper,lower,diagonal)    
    ff = -bval*u;
    bb = numeric(Nx+1);
    bb[2:(Nx)] = u[2:(Nx)]+mu*(u[1:(Nx-1)] - 2*u[2:(Nx)] + u[3:(Nx+1)]) + dt*ff[2:Nx];
    usol <- solve(A,bb)
    return(usol)
}

tvals = seq(0,T,by=dt)

uvals_FE = matrix(nrow=(Nx+1),ncol=(Nt+1))
uvals_FE[,1]=u;

truesol = matrix(nrow=(Nx+1),ncol=(Nt+1))
truesol[,1]=u;


uvals_split = matrix(nrow=(Nx+1),ncol=(Nt+1))
uvals_split[,1]=u;


for (i in 1:Nt){
    truesol[,i+1] = exact(x,tvals[i+1])
    
    uold_FE = uvals_FE[,i]
    uu_FE = diffusion_step(uold_FE,b,dt)
    uvals_FE[,i+1] = uu_FE;
    
    uold_split = uvals_split[,i]
    if (i==1){
      uuu_split = rxn_step(uold_split,dt/2)
      uu_split = diffusion_step(uuu_split,0,dt)
      u_split = rxn_step(uu_split,dt/2)
      uvals_split[,i+1] = u_split;
      }
    else {
        uprev = uvals_split[,i-1]
        uuu_split = rxn_step(uold_split,dt/2)
        uu_split = diffusion_step(uuu_split,0,dt)
        u_split = rxn_step(uu_split,dt/2)
        uvals_split[,i+1] = u_split;        
    }
}



library(plotly)
errs_FE = abs(uvals_FE-truesol)/(1+truesol)
errs_split = abs(uvals_split-truesol)/(1+truesol)


normed_err_FE = numeric(Nt+1);
normed_err_split = numeric(Nt+1);

for (i in 1:(Nt+1)){
    normed_err_FE[i] = norm(truesol[,i]-uvals_FE[,i],type="2");
    normed_err_split[i] = norm(truesol[,i]-uvals_split[,i],type="2");
}
plot(normed_err_FE,type='l')
lines(normed_err_split,col=2)



p <- plot_ly(y=tvals, x=x)%>%
    add_surface(z=t(uvals_split))
    #%add_surface(z=t(errs_split),opacity=.3)%>%
   # add_surface(z=t(errs_FE),opacity=.5)



#, z = t(uvals))%>% add_surface()
#p2 <- plot_ly(y=tvals, x=x, z = t(uvals_split))%>% add_surface()


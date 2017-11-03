# distribution, cdf, quantile and random functions for Pareto distributions
dpareto <- function(x, xm, alpha) ifelse(x > xm , alpha*xm**alpha/(x**(alpha+1)), 0)
ppareto <- function(q, xm, alpha) ifelse(q > xm , 1 - (xm/q)**alpha, 0 )
qpareto <- function(p, xm, alpha) ifelse(p < 0 | p > 1, NaN, xm*(1-p)**(-1/alpha))
rpareto <- function(n, xm, alpha) qpareto(runif(n), xm, alpha)

## n+1 cumulative sums of unif(0,1)
discrete_unif <- function(n){
  xn <- runif(n, min=0,max=1);
  yn <- c(0,cumsum(xn));
}

discrete_exp <- function(n,mu){
  xn <- rexp(n,rate=1/mu);
  yn <- c(0,cumsum(xn));
}

discrete_norm <- function(n,mu){
  xn <- rnorm(n,mean=mu,sd=1);
  yn <- c(0,cumsum(xn));
}

discrete_pareto <- function(n,xm,alpha){
  xn <- rpareto(n,xm=xm,alpha=alpha);
  yn <- c(0,cumsum(xn));
}



library(plotly)
# create data
aval <- list()
nvals = c(10, 50, 100, 250, 500, 1000)
sn <- discrete_pareto(1000,1,1.5)

l = length(nvals)

for(i in 1:l){
  aval[[i]] <-list(visible = FALSE,
                      name = paste0('N = ', nvals[i]),
                      x=0:nvals[i],
                      y=sn[1:(nvals[i]+1)])
}

aval[[1]]$visible = TRUE


# create steps and plot all traces
steps <- list()
p <- plot_ly()
for (i in 1:l) {
  p <- add_trace(p,x=aval[[i]]$x,  y=aval[[i]]$y, visible = aval[[i]]$visible, 
                 name = aval[[i]]$name, type='scatter', mode = 'markers', hoverinfo = aval[[i]]$y, 
                  showlegend = FALSE, marker=list(size=12,color='#16B1B5'))
  
  step <- list(args = list('visible', rep(FALSE, length(aval))),
               method = 'restyle',label=nvals[i])
  step$args[[2]][i] = TRUE  
  steps[[i]] = step 
}  

# add slider control to plot
p <- p %>%
  layout(sliders = list(
    list(active = 0, currentvalue = list(prefix = "N= "),
          steps = steps)
    ), xaxis=list(title="n",titlefont=list(size=12),anchor='free',position=0.0), yaxis=list(title="Sn",titlefont=list(size=12)), title = "cumulative sums Sn for Xn ~ unif(0,1)"
)


### continuous stuff


continuous_norm <- function(M,N,mu){
    sn <- discrete_norm(N,mu);
    tvals = seq(0,1,length.out=M);
    sn_cont = sn[floor(tvals*N)+1]
}


#M=1000;
#tvals <- seq(0,1,length.out=M);



aval2 <- list()
l = length(nvals)
sn <- discrete_norm(1000,0)
for(i in 1:l){
    aval2[[i]] <-list(visible = FALSE,
                     name = paste0('N = ', nvals[i]),
                     x=(0:nvals[i])/nvals[i],#tvals,
                     y=sn[1:(nvals[i]+1)])#continuous_norm(M,nvals[i],0))
}
aval2[[1]]$visible = TRUE
# create steps and plot all traces
steps2 <- list()
p2 <- plot_ly()
for (i in 1:l) {
    p2 <- add_lines(p2,x=aval2[[i]]$x,  y=aval2[[i]]$y, visible = aval2[[i]]$visible, 
                   name = 'step', line = list(shape = "hv"),
                   showlegend = TRUE,color=I('#16B1B5'))
    
    step2 <- list(args = list('visible', rep(FALSE, 3*length(aval2))),
                 method = 'restyle',label=nvals[i])
    step2$args[[2]][i] = TRUE
    step2$args[[2]][i+l] = TRUE    
    step2$args[[2]][i+2*l] = TRUE    
    steps2[[i]] = step2
}  

for (i in 1:l){
    p2 <- add_lines(p2,x=aval2[[i]]$x,  y=aval2[[i]]$y, visible = aval2[[i]]$visible, 
                    name = 'linear interp', line = list(shape = "linear",dash='dash'),
                    showlegend = TRUE,color=I('#fc8d62'))
}

for (i in 1:l){
    p2 <- add_trace(p2,x=aval2[[i]]$x,  y=aval2[[i]]$y, visible = aval2[[i]]$visible, 
                   name = 'discrete sum', type='scatter', mode = 'markers',
                   showlegend = TRUE, marker=list(size=12/(.75*i),color='#8da0cb',opacity=0.65,layer=3))
}

# add slider control to plot
p2 <- p2 %>%
    layout(sliders = list(
        list(y=-0.1,active = 0, currentvalue = list(prefix = "n= "),
             steps = steps2)), 
        xaxis=list(title="t",titlefont=list(size=12),anchor='free',position=0.0),
        yaxis=list(title="Sn(t)",titlefont=list(size=12)), 
        title = "continuous version of Sn(t) for Xn ~ unif(0,1)",
        legend = list(x = 0.15, y = 1.03,orientation = 'h', bgcolor = rgb(1,1,1,0) )
    )



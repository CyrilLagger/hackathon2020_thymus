library("ODEsensitivity")
CImod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dS <- a*S*(1-(S+R)/cMax)-q*E*S-S*(D/(D+I))-l*S-d_T*S
    dR <- a*R*(1-(S+R)/cMax)-q*E*R+l*S-d_T*R
    dE <- s+p*E*(S+R)/(g+S+R)-m*E*(S+R)-b*E*C-d_E*E
    dC <- -(k1+k2)*C+u/V1
    dD <- k12*V1/V2*C-k2*D
    return(list(c(dS, dR, dE, dC, dD)))
  })
}
CIpars <- c("a","cMax", "q","I","l","d_T","s","p","g","m","b","d_E","k1","k2","u","V1","V2","k12")
CIbinf <- c(.54,883260000,9e-8,10,10e-10,0.153,1.3e5,.11205,1.8171e7,3.0798e-10,2e-6,0.03708,1.44,.72,0,22.5,13.5,0.36)
CIbsup <- c(.66,1079540000,1.5e-7,100,10e-8,0.187,1.7e5,0.13695,2.2209e7,3.7642e-10,2e-2,0.04532,1.76,0.88,110,27.5,16.5,0.44)
CIinit <- c(S = 300, R = 0, E = 10000, C=0, D=0)
CItimes <- c(0.01, seq(1, 700, by = 1))
set.seed(59281)

CIres_sobol <- ODEsobol(mod = CImod,
                        pars = CIpars,
                        state_init = CIinit,
                        times = CItimes,
                        n = 1000,
                        rfuncs = "runif",
                        rargs = paste0("min = ", CIbinf,
                                       ", max = ", CIbsup),
                        sobol_method = "Martinez",
                        ode_method = "lsoda",
                        parallel_eval = TRUE,
                        parallel_eval_ncores = parallel::detectCores())

str(CIres_sobol, vec.len = 5, give.attr = FALSE)
plot(CIres_sobol, pars_plot = c("a","q","I","l"),state_plot = "S")
plot(CIres_sobol, pars_plot = c("a","q","l"),state_plot = "R")
plot(CIres_sobol, pars_plot = c("s","p","g","m","b"),state_plot = "E")
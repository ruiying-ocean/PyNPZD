# generate.py verifies the archive before this driver executes original code.
out <- normalizePath(commandArgs(trailingOnly=TRUE)[1])
original <- readLines("EMPOWER1.r", encoding="latin1")
evaluate <- function(first, last) eval(parse(text=original[first:last]), envir=.GlobalEnv)
# Skip the original rm(list=ls()); load its functions and initialization only.
evaluate(14, 598)
options(digits=17)
save_csv <- function(x, name) write.csv(x, file.path(out, name), row.names=FALSE)

parameter_names <- c("Vp0", "alpha", "kN", "mP", "mP2", "Imax", "kz", "phiP", "phiD", "betaz", "kNz", "mz", "mz2", "VD", "mD", "wmix", "CtoChl", "kw", "kc")
setting_names <- c("Pinit", "Ninit", "Zinit", "Dsinit", "nyears", "flag_stn", "flag_LI", "flag_atten", "flag_irrad", "flag_PIcurve", "flag_grazing")
toml_values <- function(names) paste(names, "=", sprintf("%.15g", sapply(names, get)))
writeLines(c(
  "# Generated from the official EMPOWER-1.0 inputs; independent of config.toml.",
  "[PARAMETER]", toml_values(parameter_names),
  "", "[SETUP]", toml_values(setting_names),
  paste("delta_t =", sprintf("%.15g", delta_t*24), "# hours between Python outputs"),
  "flag_force = false", "", "[FORCING]", 'source = "station"'
), file.path(out, "config.toml"))
save_csv(forcing, "station_table.csv")

forcing_cases <- list()
for (station in 1:4) {
  flag_stn <- station
  evaluate(515, 598)
  forcing_cases[[station]] <- data.frame(station=station, day=1:366, MLD=MLD, SST=SST,
    N0=aN*MLD+bN, day_length=FNdaylcalc(1:366, latradians),
    noon_par=FNnoonparcalc(1:366, latradians, clouds, e0))
}
save_csv(do.call(rbind, forcing_cases), "forcing.csv")

light <- expand.grid(depth=c(3,5,10,23,100), Iin=c(10,100,300), Vp=c(5,10), irrad=1:2, PI=1:2)
light$Iout <- light$Iin*exp(-0.08*light$depth)
light$num <- mapply(function(z,i,o,v,ir,p) FNLIcalcNum(z,i,o,0.08,0.15,v,12,ir,p),
  light$depth, light$Iin, light$Iout, light$Vp, light$irrad, light$PI)
light$ep85 <- mapply(function(z,i,o,v) FNLIcalcEP85(z,i,o,0.08,0.15,v,12),
  light$depth, light$Iin, light$Iout, light$Vp)
light$a93 <- mapply(function(z,i,o,v) FNLIcalcA93(z,i,o,0.08,0.15,v,12,0.3),
  light$depth, light$Iin, light$Iout, light$Vp)
save_csv(light, "light.csv")
absorption <- expand.grid(ss=c(0.1,0.5,1,2), top=c(0,5,23), thickness=c(3,18,100))
absorption$value <- mapply(function(s,z,h) FNaphy(s,z,z+h,0.3),
  absorption$ss, absorption$top, absorption$thickness)
save_csv(absorption, "absorption.csv")

flag_stn <- 2
evaluate(515, 598)
choices <- data.frame(LI=c(1,1,1,1,2,3), irrad=c(1,1,2,2,2,2), PI=c(1,2,1,2,1,1))
cases <- expand.grid(choice=1:6, atten=1:2, grazing=1:2, depth=c(3,5,10,23,100), change=c(-1,0,1), P=c(0.1,1.2), day=c(30,180))
cases <- cbind(choices[cases$choice,], cases[,names(cases)!="choice"])
cases$N <- 2
cases$Z <- 0.2
cases$D <- 0.3
cases$temperature <- 10
base <- which(cases$LI==1 & cases$irrad==2 & cases$PI==1 & cases$atten==2 &
  cases$grazing==2 & cases$depth==100 & cases$change==0 & cases$P==0.1 & cases$day==180)
stopifnot(length(base)==1)
edges <- cases[rep(base,8),]
edges[,c("P","N","Z","D")] <- rbind(
  c(0,0,0,0), c(0,10,0.2,0.3), c(0.1,0,0.2,0.3), c(0.1,2,0,0.3),
  c(0.1,2,0.2,0), c(1.2,0.01,0.5,1), c(0.1,2,0.2,0.3), c(0.1,2,0.2,0.3))
edges$temperature[7:8] <- c(0,25)
cases <- rbind(cases, edges)
rates <- matrix(NA_real_, nrow(cases), 4)
flux <- matrix(0,10,4)
Y <- rep(0,20)
nstepsday <- 24
istepnow <- 13
for (i in seq_len(nrow(cases))) {
  a <- cases[i,]
  flag_LI <- a$LI
  flag_irrad <- a$irrad
  flag_PIcurve <- a$PI
  flag_atten <- a$atten
  flag_grazing <- a$grazing
  daynow <- a$day
  MLD <- rep(a$depth,366)
  MLD[daynow+1] <- a$depth+a$change
  SST <- rep(a$temperature,366)
  rates[i,] <- colSums(FNget_flux(c(a$P,a$N,a$Z,a$D)))
}
colnames(rates) <- c("dP","dN","dZ","dD")
stopifnot(nrow(cases)==1448, all(is.finite(rates)))
save_csv(cbind(cases,rates), "rates.csv")

# Execute the original five-year RK4 loop, suppressing its rounded files/plots.
evaluate(461, 598)
flag_outtype <- 0
evaluate(606, 776)
stopifnot(nrow(Svar)==1826, all(is.finite(Svar)), all(Svar>=0))
save_csv(data.frame(time=tm,P=Svar[,1],N=Svar[,2],Z=Svar[,3],D=Svar[,4]), "biotrans_5years.csv")

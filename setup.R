library(cOde)
library(ggplot2)
library(tidyverse)
library(data.table)
library(deSolve)
library(dMod)
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

# add data and make useful
# data <- as.data.table(read_csv("./20220520_105706_df_mean_std.csv"))
data <- as.data.table(read_csv('./geometric_mean_nonorm.csv'))

data[,`:=`(minutes=NULL, time_point=NULL)]
value_names <- names(data)[grep("au$", names(data))]
# sigma_names <- names(data)[grep("std$", names(data))]

data_long <- melt(
  data,
  # measure.vars = list(value_names, sigma_names),
  measure.vars = list(value_names),
  # value.name = c("value", "sigma"),
  value.name = c("value"),
  variable.name = "name"
)


setnames(data_long, "seconds", "time", skip_absent = TRUE)
data_long[, condition:="same_conc"]
data_long[, name:=str_remove(value_names, "_au")[name]]
data_long[, name := paste0(name, "_obs")]

use_data <- data_long %>% as.data.frame() %>% as.datalist(split.by = "condition")

usedata_long <- as.data.table(as.data.frame(use_data))

# 2. Equation list --------------------------------------------------------

el <- NULL
el <- addReaction(el, 
                 from = "A + A", 
                 to = "AA", 
                 rate = "k_2x_1xx * A * A",
                 description = "production of AA by A + A")
el <- addReaction(el, 
                 from = "AA", 
                 to = "A + A", 
                 rate = "k_1xx_2x * AA",
                 description = "decay of AA to A + A")
el <- addReaction(el,
                 from = "A + AA",
                 to = "AAA"     ,
                 rate = "k_1x1xx_1xxx * AA * A",
                 description = "degredation of A + AA to AAA")
el <- addReaction(el,
                 from = "AAA",
                 to = "A + AA"     ,
                 rate = "k_1xxx_1x1xx * AAA",
                 description = "degredation of AAA to A + AA")

el <- addReaction(el,
                 from = "B + B",
                 to = "BB",
                 rate = "k_2x_1xx * B * B",
                 description = "production of BB by B + B")
el <- addReaction(el,
                 from = "BB",
                 to = "B + B",
                 rate = "k_1xx_2x * BB",
                 description = "decay of BB to B + B")
el <- addReaction(el,
                 from = "B + BB",
                 to = "BBB"     ,
                 rate = "k_1x1xx_1xxx * BB * B",
                 description = "degredation of B + BB to BBB")
el <- addReaction(el,
                 from = "BBB",
                 to = "B + BB"     ,
                 rate = "k_1xxx_1x1xx * BBB",
                 description = "degredation of BBB to B + BB")

el <- addReaction(el,
                 from = "A + B",
                 to = "AB"     ,
                 rate = "k_2x_1xx * A * B",
                 description = "Production of AB from A + B")

el <- addReaction(el,
                 from = "AB",
                 to = "A + B"     ,
                 rate = "k_1xx_2x * AB",
                 description = "Production of A + B from AB")
el <- addReaction(el,
                 from = "AA + B",
                 to = "AAB"     ,
                 rate = "k_1x1xx_1xxx* AA * B",
                 description = "Production of AAB from AA + B")
el <- addReaction(el,
                 from = "AAB",
                 to = "AA + B"     ,
                 rate = "k_1xxx_1x1xx * AAB",
                 description = "Production of AA + B from AAB")

el <- addReaction(el,
                 from = "A + BB",
                 to = "ABB"     ,
                 rate = "k_1x1xx_1xxx * A * BB",
                 description = "Production of ABB from A + BB")
el <- addReaction(el,
                 from = "ABB",
                 to = "A + BB"     ,
                 rate = "k_1xxx_1x1xx * ABB",
                 description = "Production of A + BB from A + BB")

el <- addReaction(el,
                 from = "AB + B",
                 to = "ABB"     ,
                 rate = "k_1x1xx_1xxx * AB * B",
                 description = "Production of ABB from AB + B")
el <- addReaction(el,
                 from = "ABB",
                 to = "AB + B"     ,
                 rate = "k_1xxx_1x1xx * ABB",
                 description = "Production of AB + B from ABB")

el <- addReaction(el,
                 from = "A + AB",
                 to = "AAB"     ,
                 rate = "k_1x1xx_1xxx * A * AB",
                 description = "Production of AAB from A + AB")
el <- addReaction(el,
                 from = "AAB",
                 to = "A + AB"     ,
                 rate = "k_1xxx_1x1xx * AAB",
                 description = "Production of A + AB from AAB")


# 3. Observables ----------------------------------------------------------
observables <- eqnvec(
  A_obs = "scale_A * A", 
  AA_obs = "scale_AA * AA", 
  AAA_obs = "scale_AAA * AAA",
  B_obs = "scale_B * B",
  BB_obs = "scale_BB * BB",
  BBB_obs = "scale_BBB * BBB",
  AB_obs = "scale_AB * AB",
  AAB_obs = "scale_AAB * AAB",
  ABB_obs = "scale_ABB * ABB"
)

# 4. Error model ----------------------------------------------------------
obsError <- names(observables)
# errors <- as.eqnvec(
#   c(#paste0("sqrt(",obsError,"**2*sigma_rel**2+sigma_abs**2)")
#     # paste0("sqrt(",obsError,"*sigma_",obsError,"_rel**2+sigma_",obsError,"_abs**2)")
#      paste0("sigma_",obsError)
#   ), names = c(obsError
#   )
# )

errors <- as.eqnvec(
     paste0("sigma_",obsError),
     names = c(obsError
  )
)

# 5. Create ODE model -----------------------------------------------------
model <- odemodel(f = el #,
                  # forcings = NULL,
                  # events = data.frame(var = "mixing",
                  #                           time = 70,
                  #                           root=NA,
                  #                           value = 1,
                  #                           method = "replace")#,
                  # # jacobian = "inz.lsodes",
                  # compile= T,
                  # fixed = NULL,
                  # modelname = "x_Pex14"
                  )

x <- Xs(
  odemodel = model,
  optionsOde = c(rtol = 1e-10, atol = 1e-10),
  optionsSens = c(rtol = 1e-10, atol = 1e-10)
  )

# 6. Observation function -------------------------------------------------
g <- Y(g = observables,
       f = el #, compile = TRUE, modelname = "g_Pex14"
       )


# 7. Compile error model --------------------------------------------------

#errPars <- getSymbols(errors)
#errPars <- errPars[which(grepl("sigma", errPars))]
e <- Y(g = errors,
      f = c(as.eqnvec(el), observables),
      states = names(observables),
      # attach.input = FALSE,
       compile = TRUE,
      # modelname = "e_Pex14"
      )

# Generate observation functions
# g <- Y(observables, x, compile = TRUE, modelname = "g_Pex14")

# 8. Transformation function ----------------------------------------------
paramInner <- c(getParameters(el), getSymbols(observables), getSymbols(errors)) %>% unique()
names(paramInner) <- paramInner

trafo <- as.eqnvec(paramInner)
conditionGrid <- attr(use_data, "condition.grid")

init_mixed <- c("AB", "AAB", "ABB")
trafo <- insert(trafo, "x~y", x="init_mixed", y=0)
trafo <- insert(trafo, "x~y", x="scale_A", y=12)
# trafo <- insert(trafo, "x~y", x="scale", y=1)
# trafo <- insert(trafo, "x~y", x="AAA", y="k_21 * k_21 / (k_12 * k_12) * A * A * A")
# trafo <- insert(trafo, "x~y", x="BBB", y="k_21 * k_21 / (k_12 * k_12) * B * B * B")
# trafo <- insert(trafo, "x~y", x="AA", y="k_21 / k_12 * A * A")
# trafo <- insert(trafo, "x~y", x="BB", y="k_21 / k_12 * B * B")


trafoList <- branch(trafo, conditionGrid)
trafoList <- insert(trafoList, "x~10^(x)", x=.currentSymbols)

p <- P(trafo = trafoList)

# 9. Compile prediction function ------------------------------------------
prdFunc <- Reduce("*", list(g, x, p))

# 10. Define Objective Function -------------------------------------------
paramOuter <- structure(rep(-1, length(getParameters(p))), names = getParameters(p))

obj <- normL2(
  data = use_data,
  x = prdFunc,
  errmodel=e
  ) + constraintL2(paramOuter, sigma = 8)

# 11. Fit -----------------------------------------------------------------
fitOutput <- mstrust(
  objfun = obj,
  center=dMod::msParframe(paramOuter, n = 100, seed=1),
  studyname = "fits",
  cores = 6,
  iterlim = 1000,
  printIter = FALSE
)

result <- as.parframe(fitOutput)

bestFit <- result[1]
10^as.parvec(bestFit)

plotValues(result)
plotCombined(prdFunc(seq(70, 700), bestFit), use_data)

plotFluxes(bestFit, prdFunc, seq(0, 700), subset(el, "A"))

# Compute the profile likelihood around the optimum
profiles <- profile(obj, bestFit, names(bestFit), limits = c(-6, 6), cores = 4)

# Take a look at each parameter
# plotProfile(profiles, mode=="data")
plotProfile(profiles, mode=="data")


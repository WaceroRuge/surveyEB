rm(list = ls())
library(pacman)
p_load(tidyverse, lme4, assertr, RColorBrewer, survey, srvyr, lmerTest, furrr, 
       progressr, ggpubr, latex2exp )
options(survey.lonely.psu = "adjust")

# Directories
inPath <- '../input'

# Data base
entra <- file.path(inPath, '0.data_dane.RData') # s = GEIH, s' = ECV
load(entra)

# Pre-process
df_geih <- df_geih |> mutate(DIRECTORIO = gsub("(.*)-(.*)-(.*)", "\\1", id))

# "Atypical" individuals
id.atipicos.geih <- c("7566491-1-2", "7555746-1-4", "7562860-1-2", "7562912-1-1",
                      "7564794-1-1", "7558151-1-1", "7554009-1-1", "7547065-1-1",
                      "7540597-1-1", "7560377-1-1", "7545920-1-1", "7547071-1-1", 
                      "7548820-1-1", "7560468-1-3", "7564696-1-2", "7564696-1-1", 
                      "7557859-1-1", "7548866-1-1", "7548619-1-1", "7540786-1-1")

id.atipicos.ecv <- NA

df_geih <- df_geih |> mutate(atipico = ifelse(id %in% id.atipicos.geih, 1, 0))
df_ecv <- df_ecv |> mutate(atipico = ifelse(id %in% id.atipicos.ecv, 1, 0))

df_geih.mod <- df_geih

# Domain sample sizes
df_geih.mod <- df_geih.mod |>
  group_by(dom) |>
  mutate(nd = n()) |>
  ungroup()

# summary(df_geih.mod$ingpcug)

# Correlation plot
cor(df_geih.mod$ingpcug, df_geih.mod$MONTO_SS)
plot(df_geih.mod$ingpcug, df_geih.mod$MONTO_SS)

# Constant for transformation
grid.const <- seq(from = 50, to = 150, by = 5)

get.skew <- function(const, df){
  library(moments)
  df_mod <- df |>
    mutate(yd = log(ingpcug + const)) %>%
    mutate(atipico.dom = ifelse(dom %in% c(27, 32), 1, 0)) |>
    mutate(MONTO_SS = log(MONTO_SS + const)) |> 
    select(yd, dom, atipico, atipico.dom, SEXO, AFILIACION_SS, MONTO_SS, NIVEL_EDU)

  mod.pob.temp <- lmer(yd ~ SEXO + AFILIACION_SS + MONTO_SS + NIVEL_EDU + 
                         atipico + atipico.dom + (1|dom), data = df_mod)

  res.temp <- residuals(mod.pob.temp)
  sk <- moments::skewness(res.temp)
  return(sk)

  ee <- resid(mod.pob.temp)
}

# sk <- map_dbl(grid.const, ~get.skew(.x, df = df_geih.mod))
# df_sk <- tibble(grid.const = grid.const, sk = sk)
# cons.skew <- df_sk[which.min(abs(df_sk$sk)),][["grid.const"]] # 730
cons.skew <- 65

doms.atip <- c(27,32)
df_geih.mod$atipico.dom <-  ifelse(df_geih.mod$dom  %in% doms.atip, 1, 0)
df_geih$atipico.dom <-  ifelse(df_geih$dom %in% doms.atip, 1, 0)
df_ecv$atipico.dom <-  ifelse(df_ecv$dom  %in% doms.atip, 1, 0)

# yd variable for model
df_geih.mod <- df_geih.mod |> 
  mutate(yd = log(ingpcug + cons.skew),
         MONTO_SS = log(MONTO_SS + cons.skew))

df_geih <- df_geih |> 
  mutate(MONTO_SS = log(MONTO_SS + cons.skew))

df_ecv <- df_ecv |> 
  mutate(MONTO_SS = log(MONTO_SS + cons.skew))

# Survey EB estimation
fit <- lmer(yd ~ SEXO + AFILIACION_SS + MONTO_SS + NIVEL_EDU + atipico + atipico.dom + (1|dom), data = df_geih.mod)

# Random effects diagnostic
uu <- ranef(fit)$dom[,1]
df.residuals <- tibble(ranef = uu)
p.qqRanef <- df.residuals |>
  ggplot(aes(sample = ranef)) +
  stat_qq() +
  stat_qq_line() +
  theme_bw() +
  ylab("Sample quantiles") +
  xlab("Theorical quantiles")

p.qqRanef

# Residual diagnostic
ee <- resid(fit)
df.residuals <- tibble(resid = ee)
p.qqRes <- df.residuals |>
  ggplot(aes(sample = resid)) +
  stat_qq() +
  stat_qq_line() +
  theme_bw() +
  ylab("Sample quantiles") +
  xlab("Theorical quantiles")

p.qqRes

p <- ggarrange(p.qqRanef, p.qqRes, ncol = 2, labels = c("A", "B"))
p

# Estimated model parameters 
betaest <- fixef(fit)
upred  <- c(ranef(fit)$dom$`(Intercept)`)
sigmau2est <- as.numeric(VarCorr(fit))
sigmae2est <- summary(fit)$sigma^2

# Direct estimators
df_geih <- df_geih |> 
  mutate(poor = ingpcug < z,
         gap = ((z-ingpcug)/z)*(ingpcug < z)) |> 
  ungroup()

diseno.geih <- svydesign(ids = ~DIRECTORIO, weights = ~fexp, strata = ~dpto,
                         data = df_geih, nest = TRUE) |> 
  as_survey()

est.DIR <- diseno.geih |> 
  group_by(dom) |> 
  summarise(nd = n(),
    est.DIR.poor = survey_mean(poor, vartype = "cv"),
    est.DIR.gap = survey_mean(gap, vartype = "cv"),
    est.DIR.income = survey_mean(ingpcug, vartype = "cv"))


# \hat{\delta_{di}} in GEIH (sample s)
df_geih <- df_geih |> 
  left_join(df_geih.mod |> select(dom, nd) |> unique(), by = "dom") |> 
  mutate(cons.skew = cons.skew,
         mudcondi = as.vector(predict(fit, newdata = df_geih)),
         gammad = sigmau2est/(sigmau2est+sigmae2est/nd)) |> 
  mutate(sigma2cond = sigmau2est*(1-gammad)+sigmae2est) |> 
  mutate(alphadi = (log(z+cons.skew)-mudcondi)/sqrt(sigma2cond))

df_geih <- df_geih |> 
  mutate(poor.delta.di = pnorm(alphadi),
         gap.delta.di = pnorm(alphadi)*(1-(1/z)*(exp(mudcondi+sigma2cond/2)*pnorm(alphadi-sqrt(sigma2cond))/pnorm(alphadi)-cons.skew )),
         income.delta.di = -cons.skew + exp(0.5 * sigma2cond + mudcondi))

diseno.geih <- svydesign(ids = ~DIRECTORIO, weights = ~fexp, strata = ~dpto,
                         data = df_geih, nest = TRUE) |> 
  as_survey()

est.SEB.geih <- diseno.geih |> 
  group_by(dom) |> 
  summarise(est.SEB_GEIH.poor = survey_mean(poor.delta.di, vartype = "cv"),
            est.SEB_GEIH.gap = survey_mean(gap.delta.di, vartype = "cv"),
            est.SEB_GEIH.income = survey_mean(income.delta.di, vartype = "cv"))

# # Prediction in s'
df_ecv$SEXO[is.na(df_ecv$SEXO)] <- sample(unique(df_ecv$SEXO), 1)
  
df_ecv <- df_ecv |>
  left_join(df_geih |> select(dom, nd, gammad) |> unique(), by = "dom") |> 
  mutate(cons.skew = cons.skew,
         mudcondi = predict(fit, newdata = df_ecv)) |> 
  mutate(sigma2cond = sigmau2est*(1-gammad)+sigmae2est) |> 
  mutate(alphadi = (log(z+cons.skew)-mudcondi)/sqrt(sigma2cond))
  
df_ecv <- df_ecv |> 
  mutate(poor.delta.di = pnorm(alphadi),
         gap.delta.di = pnorm(alphadi)*(1-(1/z)*(exp(mudcondi+sigma2cond/2)*pnorm(alphadi-sqrt(sigma2cond))/pnorm(alphadi)-cons.skew )),
         income.delta.di = -cons.skew + exp(0.5 * sigma2cond + mudcondi))


diseno.ecv <- svydesign(ids = ~MUN+DIRECTORIO, weights = ~fexp, 
          data = df_ecv,
          nest = TRUE) |> 
  as_survey()

est.SEB.ecv <- diseno.ecv |> 
  group_by(dom) |> 
  summarise(nd.prime = n(),
            est.SEB_ECV.poor = survey_mean(poor.delta.di, vartype = "cv"),
            est.SEB_ECV.gap = survey_mean(gap.delta.di, vartype = "cv"),
            est.SEB_ECV.income = survey_mean(income.delta.di, vartype = "cv"))

# Consolidate results
est.consol <- est.DIR |> 
  left_join(est.SEB.geih, by = "dom") |> 
  left_join(est.SEB.ecv, by = "dom")

# Add real domain names
id.dptos <- tribble(
  ~cod_dpto, ~dpto,
  '05', 'Antioquia','08', 'Atlántico','11', 'Bogotá, D.C.','13', 'Bolívar','15',
  'Boyacá','17', 'Caldas','18', 'Caquetá','19', 'Cauca','20', 'Cesar','23',
  'Córdoba','25', 'Cundinamarca','27', 'Chocó','41', 'Huila','44',
  'La Guajira','47', 'Magdalena','50', 'Meta','52', 'Nariño','54',
  'Norte de Santander','63', 'Quindío','66', 'Risaralda','68', 'Santander','70',
  'Sucre','73', 'Tolima','76', 'Valle del Cauca','81', 'Arauca','85',
  'Casanare','86', 'Putumayo','88', 'San Andrés','91', 'Amazonas','94',
  'Guainía','95', 'Guaviare','97', 'Vaupés','99', 'Vichada')

id.etnia <- tribble(
  ~cod_etnia, ~etnia,
  1, 'IND', 2, 'GRP', 3, 'NM', 4, 'NIN')

doms.ori.names <- expand_grid(cod_dpto = id.dptos$cod_dpto, cod_etnia = id.etnia$cod_etnia) |>
  left_join(id.dptos, by = 'cod_dpto') |> 
  left_join(id.etnia, by = 'cod_etnia') |> 
  mutate(dom.ori = paste0(cod_dpto, "-", cod_etnia),
         dom.ori.name = paste0(dpto, " - ", etnia)) |> 
  select(dom.ori, dom.ori.name)
  
df.plot <- est.consol |> 
  left_join(df_ecv |> select(dom, dom.ori) |> unique(), by = 'dom') |> 
  select(-dom) |> 
  pivot_longer(!c("dom.ori", "nd", "nd.prime")) |> 
  arrange(nd) |> 
  separate_wider_delim(name, ".", names = c("Sam", "Estimator", "Index")) |> 
  left_join(doms.ori.names, by = 'dom.ori')

temp <- df.plot |> select(dom.ori.name, nd) |> unique()
df.plot <- df.plot |> 
  mutate(dom.ori.name = factor(dom.ori.name, levels = temp$dom.ori.name, labels = temp$dom.ori.name),
         cd.prime = nd.prime / nd) |> 
  filter(!Estimator %in% c('SEB_GEIH')) |> 
  mutate(Estimator = gsub("_ECV", "", Estimator),
         value = if_else(Index %in% c("poor", "gap", "poor_cv", "gap_cv"), value*100, value))

# Estimation plot
pal <- brewer.pal(n = 11, name = "Paired")
my_theme <- theme(legend.position="top", legend.title= element_blank(), legend.direction = "horizontal",
                  axis.title.x = element_text(size = 20),
                  axis.title.y = element_text(size = 20),
                  axis.text.x = element_text(size = 13),
                  axis.text.y = element_text(size = 15),
                  legend.key.size = unit(1.8, 'cm'),
                  legend.text = element_text(size=14))


plot_estimates <- function(ind, y.lab){
  df_barras <- df.plot |> filter(Index == ind, Estimator == "DIR")  # O tomar un único set
  
  factor_escala <- max(df.plot |> filter(Index == ind) |> pull(value), na.rm = TRUE) / max(df_barras$cd.prime)
  df_barras$barras_escala <- df_barras$cd.prime * factor_escala
  
  
  
  estimates.plot <- df.plot |> 
    filter(Index == ind) |> 
    ggplot(aes(x = dom.ori.name)) +
    geom_line(aes(y = value, colour = Estimator, group = Estimator)) +
    geom_point(aes(y = value, colour = Estimator, group = Estimator, shape = Estimator), size = 2) +
    geom_bar(data = df_barras, aes(y=barras_escala), stat="identity", size=.1, 
             fill="#69b3a2", color="#69b3a2", alpha=0, group = "Estimator") +
    geom_text(data = df_barras, aes(x = dom.ori.name, y = barras_escala, label = round(cd.prime,1)), 
              vjust = 0.5, hjust = 0, size = 4, angle = 90) +
    scale_y_continuous(name = y.lab, sec.axis = sec_axis(~ ./ factor_escala,
                                          name="Sample size ratio (ECV/GEIH)")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.ticks = element_line(colour = "white"),
          axis.title.y.right = element_text(size = 14) ) +
    labs(x = "Area") +
    my_theme +
    scale_shape_manual(values=c("DIR"=16, "SEB"=17))+
    scale_color_manual(
      values =  c("DIR"=pal[1],
                  "SEB"=pal[4])) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size=6))) +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(size=20))
  
  return(estimates.plot)
  
}

plot_estimates(ind = "poor", y.lab = "Poverty rate (%)")
plot_estimates(ind = "gap", y.lab = "Poverty gap (%)")

###############################################################
# # MSE estimation
###############################################################

pb_mse_SEB_indi_Modelo_CenSin_Fast2 <- function(sam, sam2, x.vars, dom.col, sigmau2, sigmae2, cons.skew, B = 500){
  
  require(survey)
  require(srvyr)

  D <- sam |> select(dom) |> unique() |> nrow()
  
  sam2 <- sam2 |> mutate(id.temp = 1:n()) |> relocate(id.temp, .before = "dom")
  df.mse <- NULL
  
  for(b in 1:B){
    
    # Generation of bootstrap large sample:
    sam2.b <- sam2 |> select(all_of(c(dom.col, x.vars, "fexp"))) |> arrange(dom)
    sam2.b <- sam2 |> 
      group_by(dom) |> 
      mutate(ud.b  = rnorm(1, 0, sqrt(sigmau2)),
             edi.b = rnorm(nd, 0, sqrt(sigmae2))) |> 
      ungroup()
    
    Xdi.beta <- predict(fit, newdata = sam2.b, re.form = NA)
    
    
    ud.b <- sam2.b |> select(dom, ud.b) |> pull(ud.b)
    edi.b <- sam2.b |> select(dom, edi.b) |> pull(edi.b)
    sam2.b$ydi.b <- Xdi.beta + ud.b + edi.b 
    sam2.b$Edi.b <- exp(sam2.b$ydi.b) - cons.skew
    
    sam2.b <- sam2.b |> 
      mutate(deltadi.poor.b = as.numeric(Edi.b<z),
             deltadi.gap.b = (Edi.b<z)*(z-Edi.b)/z,
             deltadi.income.b = Edi.b)
    
    # Generation of bootstrap sample:
    sam.b <- sam |> select(all_of(c(dom.col, x.vars, "fexp"))) |> arrange(dom)
    sam.b <- sam.b |> 
      group_by(dom) |> 
      mutate(edi.b = rnorm(n(), 0, sqrt(sigmae2))) |> 
      left_join(sam2.b |> select(dom, ud.b) |> unique(), by = "dom")
    
    Xdi.beta <- predict(fit, newdata = sam.b, re.form = NA)
    
    
    ud.b <- sam.b |> select(dom, ud.b) |> pull(ud.b)
    edi.b <- sam.b |> select(dom, edi.b) |> pull(edi.b)
    ydi.b <- Xdi.beta + ud.b + edi.b
    sam.b$ydi.b <- ydi.b
    
    # Bootstrap model fitting and estimation
    mod.b <- lmer(ydi.b ~ SEXO + AFILIACION_SS + MONTO_SS + NIVEL_EDU + atipico + atipico.dom + (1|dom), data = sam.b)
    
    # # Empirical Best predictions
    sigmau2est.b <- as.numeric(VarCorr(mod.b))
    sigmae2est.b <- summary(mod.b)$sigma^2
    sizes.nd.sam <- sam |> select(dom, nd) |> unique()
    sizes.nd.sam2 <- sam2 |> 
      group_by(dom) |> 
      summarise(nd = n()) |> 
      select(dom, nd) |> 
      unique()
    gammad.b <- sigmau2est.b/(sigmau2est.b+sigmae2est.b/sizes.nd.sam$nd)
    
    mudcond.b <- predict(mod.b, newdata = sam2.b) |> as.vector()
    sigma2cond.b <- tibble(dom = unique(sam2.b$dom), sigma2cond.b = sigmau2est.b*(1-gammad.b) + sigmae2est.b)
    
    sam2.b$mudcond.b <- mudcond.b
    sam2.b$cons.skew <- cons.skew
    
    sam2.b <- sam2.b |> 
      left_join(sigma2cond.b, by = 'dom') |> 
      mutate(alphadi.b = (log(z+cons.skew)-mudcond.b)/sqrt(sigma2cond.b)) |> 
      mutate(deltadiEB.poor.b = pnorm(alphadi.b),
             deltadiEB.gap.b = pnorm(alphadi.b)*(1-(exp(mudcond.b+sigma2cond.b/2)*pnorm(alphadi.b-sqrt(sigma2cond.b))/pnorm(alphadi.b)-cons.skew )/z),
             deltadiEB.income.b = -cons.skew + exp(0.5 * sigma2cond.b + mudcond.b))
    
    
    # Helper function to compute estimates
    compute_estimates <- function(index, d, design) {
      formula <- as.formula(paste0("~deltadi.", index, ".b + deltadiEB.", index, ".b"))
      est <- survey::svymean(formula, design)
      vcov_est <- vcov(est)
      
      tibble(
        domain = d,
        index = index,
        deltadSEB = est[[2]], 
        delta.d.sT = est[[1]],
        var.delta.d.sT = vcov_est[1], 
        cov.SEBsT = vcov_est[2]
      )
    }
    
    # Initialize empty result
    df.Ld.b <- NULL
    
    # Loop through domains and indices
    id.doms <- sam |> select(dom) |> pull() |> unique()
    for (d in id.doms) {
      diseno.ecv <- try(sam2.b |>
        filter(dom == d) |> 
        as_survey_design(ids = c(MUN, DIRECTORIO), weights = c(fexp), 
                              nest = TRUE), silent = TRUE)
      if(class(diseno.ecv)[1] == 'try-error'){
        diseno.ecv <- sam2.b |>
                            filter(dom == d) |> 
                            as_survey_design(ids = c(DIRECTORIO), weights = c(fexp), 
                                             nest = TRUE)
      }
      
      if(d %in% c(16, 26)){
        diseno.ecv <- sam2.b |> 
          filter(dom == d) |> 
          as_survey_design(ids = c(DIRECTORIO), weights = c(fexp), strata = c(dpto),
                                 nest = TRUE) |> 
          as_survey()
      }
      
      indices <- c("poor", "gap", "income")
      df_temp <- purrr::map_dfr(indices, ~compute_estimates(.x, d, diseno.ecv))
      
      df.Ld.b <- bind_rows(df.Ld.b, df_temp)
    }
    
    # SEB estimator
    deltadSEB <- df.Ld.b$deltadSEB
    delta.d.sT <- df.Ld.b$delta.d.sT 
    
    Ad.b.Ld.b <- 2*(df.Ld.b$cov.SEBsT) - df.Ld.b$var.delta.d.sT
    
    # Bootstrap MSE estimator
    mse.est <- (deltadSEB - delta.d.sT)^2
    
    mse.est.t <- (mse.est + Ad.b.Ld.b)

    df.Ld.b <- df.Ld.b |> 
      mutate(mse.t = mse.est.t, mse.na = mse.est) |> 
      select(domain, index, mse.t, mse.na)
    
    df.mse <- df.mse |> bind_rows(df.Ld.b)
    
  }

  mse <- df.mse |> 
    group_by(domain, index) |> 
    summarise(mse.t = mean(mse.t), mse.na = mean(mse.na), .groups = "drop")
  
  return(mse)
  
}

est.mse.SEB <- pb_mse_SEB_indi_Modelo_CenSin_Fast2(
  sam = df_geih.mod,
  sam2 =  df_ecv |> select(-all_of(c("gammad", "mudcondi", "sigma2cond", "alphadi", "gap.delta.di", "poor.delta.di", "income.delta.di"))),
  x.vars = c("AFILIACION_SS", "SEXO", "MONTO_SS", "NIVEL_EDU", "atipico", "atipico.dom"),
  dom.col = "dom",
  sigmau2 = sigmau2est,
  sigmae2 = sigmae2est, 
  cons.skew = cons.skew,
  B = 500)

est.mse.SEB <- est.mse.SEB |>
  mutate(mse.cp = if_else(mse.t < 0, mse.na, mse.t),
         flag.t = if_else(mse.t < 0, 1, 0))

est.mse.SEB



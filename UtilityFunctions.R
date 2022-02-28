
# Load packages
#if (!require("groundhog")) install.packages("groundhog") #Use standard packages to ensure reproducibility
if (!require("tidyverse")) install.packages("tidyverse") #Range of convenience functions
if (!require("magrittr")) install.packages("magrittr") #Range of convenience functions
if (!require("psych")) install.packages("psych") #Range of convenience functions
if (!require("lm.beta")) install.packages("lm.beta") #Get standardized betas from regression models
if (!require("lmerTest")) install.packages("lmerTest") #Performs significance tests for LMMs
if (!require("r2glmm")) install.packages("r2glmm") #Gets effect size for LMMs
if (!require("broom")) install.packages("broom") #Helps organize analysis output
if (!require("ggeffects")) install.packages("ggeffects") #Generates predicted values from models for plotting
if (!require("scales")) install.packages("scales") #Allows rescaling of variables
if (!require("kableExtra")) install.packages("kableExtra") #For rendering tables in HTML (Rmd)
if (!require("printr")) install.packages("printr") #Efficient backup method of rendering tables in HTML
if (!require("extrafont")) install.packages("extrafont") #For expanding font options
if (!require("Cairo")) install.packages("Cairo") #For rendering graphs with high quality anti-aliasing
if (!require("sjstats")) install.packages("sjstats") #Gives CIs around Std.Beta
if (!require("Hmisc")) install.packages("Hmisc") #For corstars
if (!require("xtable")) install.packages("xtable") #For corstars
if (!require("effsize")) install.packages("effsize") #For effect sizes
if (!require("afex")) install.packages("afex") #For ANOVA
if (!require("emmeans")) install.packages("emmeans") #For ANOVA
select <- dplyr::select #Ensuring that the default for "select" is from dplyr


#### PLOTTING ####
# Set up plot theme
apatheme <- theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text = element_text(family = "Times New Roman", size = 16))

#### TABLE FUNCTION ####
freq <- function(x, digits = 2, useNA = "ifany"){
  table <- as.data.frame(table(x, useNA = useNA))
  prop <- as.data.frame(prop.table(table(x, useNA = useNA))) %>% 
    mutate(`%` = round(100*Freq,digits)) %>% select(-Freq)
  final <- left_join(table, prop, by="x") %>% dplyr::rename(Response = x)
  return(final)
}

#### REGRESSION FUNCTIONS ####
display.lm <- function(model, digits = 2, 
                       std = FALSE, CI = TRUE, conf.level = .95) {
  df <- summary(model) %>% with(df)
  display.obj <- tidy(model, conf.int = CI) %>% 
    as.data.frame() %>% mutate(p.value.char = ifelse(p.value < .001, "< .001", p.value),
                               p.value.3 = round(p.value, 3),
                               ` ` = ifelse(p.value < .001, "***", 
                                            ifelse(p.value < .01, "**", 
                                                   ifelse(p.value < .05, "\\*", 
                                                          ifelse(p.value < .10, "\\.", " "))))) %>%
    mutate(p.value = as.character(ifelse(p.value < .001, "< .001", p.value.3)),
           df = df[2])
  if(CI == T){
    display.obj <- display.obj %>% 
      dplyr::select(term, B = estimate, SE = std.error, t = statistic, df, 
                    p.value, conf.low, conf.high, ` `)
  } else {
    display.obj <- display.obj %>% 
      dplyr::select(term, B = estimate, SE = std.error, t = statistic, df, 
                    p.value, ` `)
  }
  if(std == T){
    display.obj <- display.obj %>% 
      mutate(Std.Beta = lm.beta(model) %>% coefficients() %>% as.vector()) %>% 
      dplyr::select(term, B, SE, Std.Beta, everything())
  } else {display.obj <- display.obj}
  display.obj <- display.obj %>% 
    mutate_if(is.numeric, round, digits)
  return(display.obj)
}


## Output results of a linear model using ANOVA statistics
lm_to_anova <- function(model, digits = 2, type = "III") {
  
  anova.out <- Anova(model, type = type) 
  es <- eta_squared(anova.out) %>% as.data.frame() %>% 
    select(Effect = Parameter, pes = Eta2_partial)
  prelim.df <- as.data.frame(anova.out) %>% rownames_to_column("Effect")
  resid <- prelim.df %>% filter(Effect == "Residuals") %>% 
    mutate(MSE = `Sum Sq` / Df)
  display.obj <- prelim.df %>% 
    mutate(df = str_c(Df, ", ", resid$Df),
           MSE = resid$MSE, 
           p.value = ifelse(`Pr(>F)` < .001, "< .001", as.character(round(`Pr(>F)`, 3))),
           ` ` = ifelse(`Pr(>F)` < .001, "***", 
                        ifelse(`Pr(>F)` < .01, "**", 
                               ifelse(`Pr(>F)` < .05, "\\*", 
                                      ifelse(`Pr(>F)` < .10, "\\.", " "))))) %>% 
    filter(Effect != "Residuals" & Effect != "(Intercept)") %>% 
    left_join(es) %>% 
    mutate_if(is.numeric, round, digits) %>% 
    select(Effect, df, MSE, `F` = `F value`, pes, p.value, ` `)
  return(display.obj)
}

display.lmer <- function(model, digits = 2, CI = TRUE) {
  display.obj <- summary(model) %>% coefficients() %>% as.data.frame() %>% 
    rownames_to_column() %>% rename(term = rowname, p.value = `Pr(>|t|)`) %>% 
    as.data.frame() %>% mutate(p.value.char = ifelse(p.value < .001, "< .001", p.value),
                               p.value.3 = round(p.value, 3),
                               ` ` = ifelse(p.value < .001, "***", 
                                            ifelse(p.value < .01, "**", 
                                                   ifelse(p.value < .05, "\\*", 
                                                          ifelse(p.value < .10, "\\.", " "))))) %>%
    mutate(p.value = as.character(ifelse(p.value < .001, "< .001", p.value.3))) %>%
    select(term, Estimate, SE = `Std. Error`, t = `t value`, df, p.value, ` `) %>%
    mutate_if(is.numeric, round, digits)
  if(CI == TRUE) {
    CI <- confint(model) %>% as.data.frame() %>% rownames_to_column("term") %>% 
      rename(conf.low = `2.5 %`, conf.high = `97.5 %`) %>% 
      mutate(conf.low = round(conf.low, digits), conf.high = round(conf.high, digits))
    display.obj.print <- left_join(display.obj, CI) %>% 
      select(term, Estimate, SE, t, df, p.value, conf.low, conf.high, ` `)
  } else {
    display.obj.print <- display.obj
  }
  return(display.obj.print)
}


random.lmer <- function(model, digits = 2) {
  display.obj <- VarCorr(model) %>% as.data.frame() %>% 
    select(Groups = grp, Var1 = var1, Var2 = var2, 
           `Variance/\nCovariance` = vcov, 
           `SD/Correlation` = sdcor) %>%
    mutate_if(is.numeric, round, digits) %>% 
    mutate(Var1 = ifelse(is.na(Var1) == T, " ", Var1),
           Var2 = ifelse(is.na(Var2) == T, " ", Var2))
  return(display.obj)
}


display.glm <- function(model, digits = 2, 
                        OR = FALSE, CI = TRUE, conf.level = .95) {
  df <- summary(model) %>% with(df)
  display.obj <- tidy(model, conf.int = CI) %>% 
    as.data.frame() %>% mutate(p.value.char = ifelse(p.value < .001, "< .001", p.value),
                               p.value.3 = round(p.value, 3),
                               ` ` = ifelse(p.value < .001, "***", 
                                            ifelse(p.value < .01, "**", 
                                                   ifelse(p.value < .05, "\\*", 
                                                          ifelse(p.value < .10, "\\.", " "))))) %>%
    mutate(p.value = as.character(ifelse(p.value < .001, "< .001", p.value.3)),
           df = df[2])
  if(OR == T){
    if(CI == T){
      odds <- exp(cbind(OR = coef(model), confint(model, level = conf.level))) %>% 
        as.data.frame() %>% rownames_to_column("term")
      names(odds) <- c("term", "OR", "conf.low", "conf.high")
      
      display.obj <- display.obj %>% left_join(odds) %>% 
        dplyr::select(term, B = estimate, SE = std.error, Odds.Ratio = OR, z = statistic, df, 
                      p.value, conf.low, conf.high, ` `)
    } else {
      odds <- exp(coef(model)) %>% as.vector()
      display.obj <- display.obj %>% 
        mutate(OR = exp(coef(model)) %>% as.vector()) %>% 
        dplyr::select(term, B = estimate, SE = std.error, Odds.Ratio = OR,
                      z = statistic, df, p.value, ` `)
    }
  } else {
    if(CI == T){
      display.obj <- display.obj %>% 
        dplyr::select(term, B = estimate, SE = std.error, z = statistic, df, 
                      p.value, conf.low, conf.high, ` `)
    } else {
      display.obj <- display.obj %>% 
        dplyr::select(term, B = estimate, SE = std.error, z = statistic, df, 
                      p.value, ` `)
    }
  }
  display.obj <- display.obj %>% 
    mutate_if(is.numeric, round, digits)
  return(display.obj)
}


fixed.glmer <- function(model, digits = 2) {
  df <- summary(model) %>% coefficients() %>% as.data.frame() %>% rownames_to_column()
  CI <- confint(model, parm="beta_", method="Wald") %>% 
    as.data.frame() %>% rownames_to_column() %>% 
    rename(conf.low = `2.5 %`, conf.high = `97.5 %`)
  display.obj <- left_join(df, CI) %>% rename(term = rowname, p.value = `Pr(>|z|)`) %>% 
    as.data.frame() %>% mutate(p.value.char = ifelse(p.value < .001, "< .001", p.value),
                               p.value.3 = round(p.value, 3),
                               ` ` = ifelse(p.value < .001, "***", 
                                            ifelse(p.value < .01, "**", 
                                                   ifelse(p.value < .05, "\\*", 
                                                          ifelse(p.value < .10, "\\.", " "))))) %>%
    mutate(p.value = as.character(ifelse(p.value < .001, "< .001", p.value.3))) %>%
    select(term, Estimate, SE = `Std. Error`, Z = `z value`, p.value, conf.low, conf.high, ` `) %>%
    mutate_if(is.numeric, round, digits)
  return(display.obj)
}



## Additional steps for formatting kable
htable <- function(display.obj, style = "default") {
  if (style == "default"){kab <- kable(display.obj)}
  if (style == "center"){kab <- kable(display.obj, align = c("l", rep('c', ncol(display.obj) - 1)))}
  kab %>% kable_styling(bootstrap_options = c("hover", "condensed", "responsive"), 
                        full_width = F, position = "center")
}

## SIMPLE EFFECTS
simple.effects <- function(model, x, m, minmax, by = 1) {
  terms <- names(model$model) #Extract names of variables in model
  df <- get(as.character(model$call)[3]) #Load dataframe used in model
  simple.list <- list() #Set up empty list for saving iterative model results later
  
  if(missing(x)) {x <- terms[2]} #Assigning first predictor of model as X if unspecified
  if(missing(m)) {m <- terms[3]} #Assigning second predictor in model as M if unspecified
  
  #Accounting for potentially dichotomous X variable (dummy codes it)
  if(!is.numeric(pull(df[x]))) {df[x] <- as.numeric(factor(pull(df[x]))) - 1}
  
  #For continuous moderators (i.e., > 2 levels)
  level.check <- pull(df[m]) %>% na.omit()
  if(!is.character(level.check) & length(unique(level.check)) > 2) { 
    m.vector <- pull(df[m]) #Isolating moderator as vector
    
    #Identifying levels of m at which to test effect of x
    if(missing(minmax)) {minmax <- c(floor(min(df[m], na.rm = T)), 
                                     ceiling(max(df[m], na.rm = T)))} #Assigning minmax from m's min and max values, rounded down and up to the nearest whole number, respectively
    if(is.character(minmax) == T) {
      if(minmax == "meansd") {
        #If user chose meansd, use -1SD, Mean, +1SD as levels of m for simple effects
        levels <- c(mean(m.vector, na.rm = T)-sd(m.vector, na.rm = T),
                    mean(m.vector, na.rm = T),
                    mean(m.vector, na.rm = T)+sd(m.vector, na.rm = T)) 
      } else {
        levels <- seq(floor(min(df[m], na.rm = T)), ceiling(max(df[m], na.rm = T)), by)
        print("NOTE: Unrecognized minmax type. Using default levels of moderator.")
      }
    } else {
      if (length(minmax) == 2) {
        levels <- seq(minmax[1], minmax[2], by)
      } else {
        levels <- minmax
      }
    }
    
    #Save all levels of m alongside whole numbers (i) for use in the list
    list.levels <- data.frame(levels) %>% rownames_to_column("i") %>% 
      mutate(i = as.numeric(i))
    
    #Run a model for each level of m and extract x's effect
    for (i in list.levels$i) {
      df$m.simple <- m.vector - list.levels$levels[i]
      new.formula.chr <- str_replace(formula(model), paste("\\* ", m, sep = ""), "* m.simple")
      model.simple <- lm(as.formula(paste(new.formula.chr[2], 
                                          new.formula.chr[1], 
                                          new.formula.chr[3])), df)
      tidy.model <- display.lm(model.simple)
      simple.list[[i]] <- filter(tidy.model, term == x) %>% 
        mutate(level = list.levels$levels[i])
      
      # Consolidate results
      simple <- do.call(rbind, simple.list) %>% 
        select(level, everything())
      
      # Check for levels of moderator outside range
      check <- simple$level >= min(m.vector, na.rm = T) & 
        simple$level <= max(m.vector, na.rm = T)
      if(any(check == F)) {
        print("NOTE: One or more levels of moderator fall outside range of that variable in the data")
        simple <- mutate(simple, outside.range = check == F)
      }
    }
  } else {
    levels <- as.character(unique(level.check))
    list.levels <- data.frame(levels) %>% rownames_to_column("i") %>% 
      mutate(levels = as.character(levels), i = as.numeric(i))
    for(i in list.levels$i) {
      df <- mutate(df, m.simple = 1)
      df$m.simple[pull(df[m]) == list.levels$levels[i]] <- 0      
      new.formula.chr <- str_replace(formula(model), paste("\\* ", m, sep = ""), "* m.simple")
      model.simple <- lm(as.formula(paste(new.formula.chr[2], 
                                          new.formula.chr[1], 
                                          new.formula.chr[3])), df)
      tidy.model <- display.lm(model.simple)
      simple.list[[i]] <- filter(tidy.model, term == x) %>% 
        mutate(level = as.character(list.levels$levels[i]))
    }
    
    # Consolidate results
    simple <- do.call(rbind, simple.list) %>% 
      select(level, everything())
  }
  
  names(simple) <- str_replace(names(simple), "level", m)
  return(simple)
}

simple.2ways <- function(model, x, m, w, minmax, by = 1) {
  terms <- names(model$model) #Extract names of variables in model
  df <- get(as.character(model$call)[3]) #Load dataframe used in model
  simple.list <- list() #Set up empty list for saving iterative model results later
  
  if(missing(x)) {x <- terms[2]} #Assigning first predictor of model as X if unspecified
  if(missing(w)) {w <- terms[3]} #Assigning second predictor in model as W if unspecified
  if(missing(m)) {m <- terms[4]} #Assigning third predictor in model as M if unspecified
  
  #Accounting for potentially dichotomous X and W variables
  if(!is.numeric(pull(df[x]))) {df[x] <- as.numeric(factor(pull(df[x]))) - 1}
  if(!is.numeric(pull(df[w]))) {df[w] <- as.numeric(factor(pull(df[w]))) - 1}
  
  #For continuous moderators (i.e., > 2 levels)
  level.check <- pull(df[m]) %>% na.omit()
  if(!is.character(level.check) & length(unique(level.check)) > 2) { 
    m.vector <- pull(df[m]) #Isolating moderator as vector
    
    #Identifying levels of m at which to test effect of x
    if(missing(minmax)) {minmax <- c(floor(min(df[m], na.rm  = T)), 
                                     ceiling(max(df[m], na.rm  = T)))} #Assigning minmax from m's min and max values, rounded down and up to the nearest whole number, respectively
    if(is.character(minmax) == T) {
      if(minmax == "meansd") {
        #If user chose meansd, use -1SD, Mean, +1SD as levels of m for simple effects
        levels <- c(mean(m.vector, na.rm = T)-sd(m.vector, na.rm = T),
                    mean(m.vector, na.rm = T),
                    mean(m.vector, na.rm = T)+sd(m.vector, na.rm = T)) 
      } else {
        levels <- seq(floor(min(df[m], na.rm  = T)), 
                      ceiling(max(df[m], na.rm  = T)), by)
        print("NOTE: Unrecognized minmax type. Using default levels of moderator.")
      }
    } else {
      if (length(minmax) == 2) {
        levels <- seq(minmax[1], minmax[2], by)
      } else {
        levels <- minmax
      }
    }
    
    #Save all levels of m alongside whole numbers (i) for use in the list
    list.levels <- data.frame(levels) %>% rownames_to_column("i") %>% 
      mutate(i = as.numeric(i))
    
    #Run a model for each level of m and extract x's effect
    for (i in list.levels$i) {
      df$m.simple <- m.vector - list.levels$levels[i]
      new.formula.chr <- str_replace(formula(model), paste("\\* ", m, sep = ""), "* m.simple")
      model.simple <- lm(as.formula(paste(new.formula.chr[2], 
                                          new.formula.chr[1], 
                                          new.formula.chr[3])), df)
      tidy.model <- display.lm(model.simple)
      simple.list[[i]] <- filter(tidy.model, term == paste(x,":",w,sep="")) %>% 
        mutate(level = list.levels$levels[i])
      
      # Consolidate results
      simple <- do.call(rbind, simple.list) %>% 
        select(level, everything())
      
      # Check for levels of moderator outside range
      check <- simple$level >= min(m.vector, na.rm = T) & 
        simple$level <= max(m.vector, na.rm = T)
      if(any(check == F)) {
        print("NOTE: One or more levels of moderator fall outside range of that variable in the data")
        simple <- mutate(simple, outside.range = check == F)
      }
    }
  } else {
    levels <- as.character(unique(level.check))
    list.levels <- data.frame(levels) %>% rownames_to_column("i") %>% 
      mutate(levels = as.character(levels), i = as.numeric(i))
    for(i in list.levels$i) {
      df <- mutate(df, m.simple = 1)
      df$m.simple[pull(df[m]) == list.levels$levels[i]] <- 0      
      new.formula.chr <- str_replace(formula(model), paste("\\* ", m, sep = ""), "* m.simple")
      model.simple <- lm(as.formula(paste(new.formula.chr[2], 
                                          new.formula.chr[1], 
                                          new.formula.chr[3])), df)
      tidy.model <- display.lm(model.simple)
      simple.list[[i]] <- filter(tidy.model, term == paste(x,":",w,sep="")) %>% 
        mutate(level = as.character(list.levels$levels[i]))
    }
    
    # Consolidate results
    simple <- do.call(rbind, simple.list) %>% 
      select(level, everything())
  }
  
  names(simple) <- str_replace(names(simple), "level", m)
  return(simple)
}

#### CORRELATION FUNCTIONS ####
corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower"),
                    result=c("none", "html", "latex"), round = 2){
  #Compute correlation matrix
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value 
  
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", 
                    ifelse(p < .01, "** ", 
                           ifelse(p < .05, "*  ", 
                                  ifelse(p < .10, "+  ", "   "))))
  
  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), round))[,-1]
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  
    ## remove last column and return the correlation matrix
    Rnew <- cbind(Rnew[1:length(Rnew)-1])
    if (result[1]=="none") return(Rnew)
    else{
      if(result[1]=="html") print(xtable(Rnew), type="html")
      else print(xtable(Rnew), type="latex") 
    }
  }
  
  ## remove lower triangle of correlation matrix
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)

    ## remove last row and return the correlation matrix
    #Rnew <- cbind(Rnew[1:length(Rnew)-1])
    if (result[1]=="none") return(Rnew)
    else{
      if(result[1]=="html") print(xtable(Rnew), type="html")
      else print(xtable(Rnew), type="latex") 
    }
  }
  

} 

display.r <- function(model, digits = 2) {
  broom::tidy(model) %>% 
    select(r = estimate, df = parameter, t = statistic, p.value, conf.low, conf.high) %>% 
    mutate(p.value = ifelse(p.value < .001, "< .001", as.character(round(p.value)))) %>% 
    mutate_if(is.numeric, round, digits)
}

## T.Tests
display.t <- function(model, digits = 2, CI = FALSE){
  if(str_detect(model$method, "Two Sample t-test") == TRUE) {
    data <- unlist(str_split(model$data.name, " by ")) %>% 
      str_split("\\$")
    DV <- pull(get(data[[1]][1])[data[[1]][2]])
    Group <- pull(get(data[[2]][1])[data[[2]][2]])
    d <- effsize::cohen.d(DV, Group)
    desc.out <- psych::describeBy(DV, group = Group, mat = T) %>% 
      select(group1, M=mean, SD=sd) %>% gather(var, value, -group1) %>% 
      arrange(group1) %>% 
      mutate(Stat = str_c(var, " (", group1,")")) %>% select(-group1, -var) %>% 
      mutate(Stat = factor(Stat, levels = .$Stat[1:4])) %>% 
      spread(Stat, value)
  } else {
    if(str_detect(model$method, "Paired t-test") == TRUE) {
      
      data <- unlist(str_split(model$data.name, " and ")) %>% 
        str_split("\\$")
      df.tmp <- data.frame(Var1 = pull(get(data[[1]][1])[data[[1]][2]]), 
                           Var2 = pull(get(data[[2]][1])[data[[2]][2]]))
      names(df.tmp) <- c(data[[1]][2], data[[2]][2])
      d <- effsize::cohen.d(df.tmp[,1], df.tmp[,2], paired = TRUE, na.rm = TRUE)
      
      desc.out <-  psych::describe(df.tmp) %>% as.data.frame() %>% rownames_to_column("group1") %>% 
        select(group1, M=mean, SD=sd) %>% gather(var, value, -group1) %>% 
        arrange(group1) %>% 
        mutate(Stat = str_c(var, " (", group1,")")) %>% select(-group1, -var) %>% 
        mutate(Stat = factor(Stat, levels = .$Stat[1:4])) %>% 
        spread(Stat, value)
    } else {
      stop('The t-test method is not supported.')
    }
  }
  
  model.out <- broom::tidy(model) %>% 
    mutate(`Cohen's d` = d$estimate,
           p.value = ifelse(p.value < .001, "< .001", round(p.value, digits))) %>% 
    select(t = statistic, df = parameter, p.value, `Cohen's d`)
  
  if(CI == TRUE){
    display.obj <- bind_cols(desc.out, model.out) %>% 
      mutate(conf.low = as.vector(d$conf.int[1]),
             conf.high = as.vector(d$conf.int[2])) %>% 
      mutate_if(is.numeric, round, digits)
  } else {
    display.obj <- bind_cols(desc.out, model.out) %>% 
      mutate_if(is.numeric, round, digits)
  }
  
  return(display.obj)
}


#Alpha
display.alpha <- function(df, digits = 3, check.keys = T){
  display <- as.data.frame(df) %>% 
    psych::alpha(check.keys = check.keys) %>% with(total) %>% round(digits) %>% 
    dplyr::select(raw.alpha = raw_alpha, std.alpha, mean, sd)
  return(display)
}

#### READING AND DISPLAYING DATA ####

readQ <- function(x) {
  df.tmp <- read_csv(x, skip = 3, col_names = names(read_csv(x)))
  return(df.tmp)
}

itemsQ <- function(x) {
  items.full <- read_csv(x)
  items <- items.full[1,] %>% gather(item, text)
  return(items)
}

display.items <- function(x, item.names, simplify = T, simplify.n = 5) {
  item.list <- list()
  for (i in 1:length(item.names)) {
    item.list[[i]] <- filter(x, str_detect(item, item.names[i]))
  }
  item.out <- do.call("rbind", item.list)
  
  if(simplify == T){
    items.split <- as.data.frame(str_split(item.out$text, " ", simplify = T)) %>% 
      mutate_all(as.character)
    keep <- items.split%>%
      summarise_all(~n_distinct(.)) %>%
      select_if(. != 1) %>% 
      colnames()
    keep.df <- select(items.split, all_of(keep))
    if(ncol(items.split) - ncol(keep.df) >= simplify.n){
      final.items <- keep.df } else {
        final.items <- items.split
      }
    
    final.items[final.items==""] <- NA
    final.items <- final.items %>% 
      unite("combo", na.rm = TRUE, remove = FALSE, sep=" ")
    item.out.final <- item.out %>% mutate(text = final.items$combo)
  } else {
    item.out.final <- item.out
  }
  return(item.out.final)
}


### Allowing describeBy to accept three IVs
describeBy3 <- function(y, group) {
  group[[1]] <- as.character(group[[1]])
  group[[2]] <- as.character(group[[2]])
  group[[3]] <- as.character(group[[3]])
  groups <- do.call(rbind, group) %>% t() %>% as.data.frame() %>% 
    mutate(groupZZZ = str_c(V1, V2, V3, sep = "%"))
  
  out <- psych::describeBy(y, groups$groupZZZ, mat = T) %>% 
    as.data.frame() %>% 
    separate(group1, sep = "%", 
             into = c("group1", "group2", "group3"))
  return(out)
}


  

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
if (!require("effectsize")) install.packages("effectsize") #For effect sizes
if (!require("afex")) install.packages("afex") #For ANOVA
if (!require("emmeans")) install.packages("emmeans") #For ANOVA
if (!require("pwr")) install.packages("pwr") #For power analysis
select <- dplyr::select #Ensuring that the default for "select" is from dplyr


### Round p-values for reporting
p.round <- function(p.value, digits = 2, graded = TRUE) {
  if(digits < 2){
    print("Warning: Rounding to less than 2 digits is not recommended.")
  }
  if(digits == 2) {
    if(graded == TRUE){
      rounded.p.value <- ifelse(p.value < .001, "< .001",
             ifelse(p.value < .01,
                    sprintf("%.3f", p.value),
                    sprintf("%.2f", p.value)))
      return(rounded.p.value)
    } else {
      sprint.call <- paste0("%.", digits, "f")
      less.than <- as.numeric(paste0(".", paste0(c(rep(0,digits-1), 1), collapse = "")))
      rounded.p.value <- ifelse(p.value < less.than, 
             paste0("< ", less.than),
             sprintf(sprint.call, p.value))
      return(rounded.p.value)
    }
  }
  if(digits > 2){
    sprint.call <- paste0("%.", digits, "f")
    less.than <- as.numeric(paste0(".", paste0(c(rep(0,digits-1), 1), collapse = "")))
    rounded.p.value <- ifelse(p.value < less.than, 
           paste0("< ", less.than),
           sprintf(sprint.call, p.value))
    return(rounded.p.value)
  }
}


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
display.lm <- function(model, digits = 2, graded = TRUE, 
                       std = FALSE, CI = TRUE, conf.level = .95) {
  df <- summary(model) %>% with(df)
  display.obj <- tidy(model, conf.int = CI) %>% 
    as.data.frame() %>% mutate(p.value.3 = round(p.value, 3),
                               ` ` = ifelse(p.value < .001, "***", 
                                            ifelse(p.value < .01, "**", 
                                                   ifelse(p.value < .05, "\\*", 
                                                          ifelse(p.value < .10, "\\.", " "))))) %>%
    mutate(p.value = p.round(p.value, digits, graded),
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
lm_to_anova <- function(model, digits = 2, graded = TRUE, type = "III") {
  
  anova.out <- Anova(model, type = type) 
  es <- eta_squared(anova.out) %>% as.data.frame() %>% 
    select(Effect = Parameter, pes = Eta2_partial)
  prelim.df <- as.data.frame(anova.out) %>% rownames_to_column("Effect")
  resid <- prelim.df %>% filter(Effect == "Residuals") %>% 
    mutate(MSE = `Sum Sq` / Df)
  display.obj <- prelim.df %>% 
    mutate(df = str_c(Df, ", ", resid$Df),
           MSE = resid$MSE, 
           p.value = p.round(`Pr(>F)`, digits, graded),
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

display.lmer <- function(model, digits = 2, graded = TRUE, CI = TRUE) {
  display.obj <- summary(model) %>% coefficients() %>% as.data.frame() %>% 
    rownames_to_column() %>% rename(term = rowname, p.value = `Pr(>|t|)`) %>% 
    as.data.frame() %>% mutate(p.value.3 = round(p.value, 3),
                               ` ` = ifelse(p.value < .001, "***", 
                                            ifelse(p.value < .01, "**", 
                                                   ifelse(p.value < .05, "\\*", 
                                                          ifelse(p.value < .10, "\\.", " "))))) %>%
    mutate(p.value = p.round(p.value, digits, graded)) %>%
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


display.glm <- function(model, digits = 2, graded = TRUE, 
                        OR = FALSE, CI = TRUE, conf.level = .95) {
  df <- summary(model) %>% with(df)
  display.obj <- tidy(model, conf.int = CI) %>% 
    as.data.frame() %>% mutate(p.value.3 = round(p.value, 3),
                               ` ` = ifelse(p.value < .001, "***", 
                                            ifelse(p.value < .01, "**", 
                                                   ifelse(p.value < .05, "\\*", 
                                                          ifelse(p.value < .10, "\\.", " "))))) %>%
    mutate(p.value = p.round(p.value, digits, graded),
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


display.glmer <- function(model, digits = 2, graded = TRUE) {
  df <- summary(model) %>% coefficients() %>% as.data.frame() %>% rownames_to_column()
  CI <- confint(model, parm="beta_", method="Wald") %>% 
    as.data.frame() %>% rownames_to_column() %>% 
    rename(conf.low = `2.5 %`, conf.high = `97.5 %`)
  display.obj <- left_join(df, CI) %>% rename(term = rowname, p.value = `Pr(>|z|)`) %>% 
    as.data.frame() %>% mutate(p.value.3 = round(p.value, 3),
                               ` ` = ifelse(p.value < .001, "***", 
                                            ifelse(p.value < .01, "**", 
                                                   ifelse(p.value < .05, "\\*", 
                                                          ifelse(p.value < .10, "\\.", " "))))) %>%
    mutate(p.value = p.round(p.value, digits, graded)) %>%
    select(term, Estimate, SE = `Std. Error`, Z = `z value`, p.value, conf.low, conf.high, ` `) %>%
    mutate_if(is.numeric, round, digits)
  return(display.obj)
}


## Additional steps for formatting kable
htable <- function(display.obj, style = "default", caption = NA) {
  if (style == "default"){
    if(is.na(caption)) {kab <- kable(display.obj)}
    else {kab <- kable(display.obj, caption = caption)}}
  if (style == "center"){
    if(is.na(caption)) {kab <- kable(display.obj, align = c("l", rep('c', ncol(display.obj) - 1)))}
    else {kab <- kable(display.obj, align = c("l", rep('c', ncol(display.obj) - 1)), caption = caption)}}
  kab %>% kable_styling(bootstrap_options = c("hover", "condensed", "responsive"), 
                        full_width = F, position = "center")
}


#### CORRELATION FUNCTIONS ####
corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower"),
                    result=c("none", "html", "latex"), digits = 2){
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
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), digits))[,-1]
  
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

display.r <- function(model, digits = 2, graded = TRUE) {
  broom::tidy(model) %>% 
    select(r = estimate, df = parameter, t = statistic, p.value, conf.low, conf.high) %>% 
    mutate(p.value = p.round(p.value, digits, graded)) %>% 
    mutate_if(is.numeric, round, digits)
}

## T.Tests
display.t <- function(model, digits = 2, graded = TRUE, CI = FALSE){
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
           p.value = p.round(p.value, digits, graded)) %>% 
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


  
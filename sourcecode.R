library(tidyverse)
library(magrittr)
library(DataCombine)
library(car)
library(httr)

data <- read_csv("https://raw.githubusercontent.com/Giant316/bayesianinference_HIVmortality/master/input_files/HIVdf.csv")
data <- read_csv("HIVdf.csv")
# Compute 4 inequality measures by the definition in chap 2.2.3
Theil <- data %$% left_join(., setNames(aggregate(PerCapInc*Pop, 
by = list(CountyName = CountyName), sum), c("CountyName","IncCounty")), # aggregate the per-cap income by County
by = c("CountyName")) %>% mutate(IncZCTA = PerCapInc*Pop, zctaCount = 1, # create zctaCount for T4
IncZCTA.CPI = PerCapInc*Pop/CPI2000, cInc=IncZCTA/IncCounty, # Inc.CPI-adjusted/ZCTA & Inc/ZCTA
IncNY = sum(IncZCTA), IncNY.CPI = sum(IncZCTA.CPI), sInc= IncZCTA/IncNY, # total income share in NY, ZCTA income in relation to NY income
sInc.CPI = IncZCTA.CPI/IncNY.CPI) %$%  # CPI-adjusted sINC
left_join(., setNames(aggregate(Pop, by = list(CountyName = CountyName), sum), # agrregate population by county
c("CountyName","PopCounty")), by = c("CountyName")) %>% 
mutate(cWe = Pop/PopCounty, sWe = Pop/sum(Pop),
t1j = sInc*log(sInc/sWe), t2j = cInc*log(cInc/cWe),                  #calculate T1 & T2
t3j = sInc.CPI*log(sInc.CPI/sWe), t4j = ifelse(t2j < 0,1,0)) %$%     #calculate T3 & partial T4
left_join(setNames(aggregate(zctaCount, by = list(CountyName), sum), c("CountyName", "zctaCount")),
left_join(setNames(aggregate(t3j, by = list(CountyName), sum), c("CountyName", "T3")), # nested join operation for all Theil indexes
left_join(setNames(aggregate(t4j, by = list(CountyName), sum), c("CountyName", "T4")),
left_join(setNames(aggregate(t1j, by = list(CountyName),sum), c("CountyName", "T1")),
setNames(aggregate(t2j, by = list(CountyName),sum), c("CountyName", "T2")), 
by = c("CountyName")), by = c("CountyName")), by = c("CountyName")), 
by = c("CountyName")) %>%       #nested join of T1+T2 -> T4 -> T3
InsertRow(., c("Hamilton",1,0,0,0,0),21) %>% mutate_at(c(-1), as.numeric) %>% #insert Hamilton row, (insert row function coerce cols to char type); convert cols back to numeric
mutate(T4 = T4/zctaCount) %>% select(CountyName, T1,T2,T3,T4) # re-order sequence

# create T1 Boxplot 
Boxplot(Theil$T1, id = list(labels = Theil$CountyName, cex = 3, n = 2), 
        ann = FALSE, cex.axis = 2.5, cex = 3)
title(ylab = parse(text=paste("T", "^1 ","*Index")), cex.lab = 3, line = 4)

# create T2 Boxplot 
Boxplot(Theil$T2, id = list(labels = Theil$CountyName, cex = 3), 
        ann = FALSE, cex.axis = 2.5, cex = 3)
title(ylab = parse(text = paste("T^2 ", "*Index")), cex.lab = 3, line = 4)

# create T3 Boxplot 
Boxplot(Theil$T3, id = list(labels = Theil$CountyName, cex = 3, n = 2), 
        ann = FALSE, cex.axis = 2.5, cex = 3)
title(ylab = parse(text = paste("T^3 *Index")), cex.lab = 3, line = 4)

# create T4 Boxplot 
Boxplot(Theil$T4, id = list(labels = Theil$CountyName, cex = 3), 
        ann = FALSE, cex.axis = 2.5, cex = 3)
title(ylab = parse(text = paste("T^4 *Index")), cex.lab = 3, line = 4)

# Retrieve Poverty Info 
pov <- data %>% select(CountyName, PovertyPercentage) %>% unique() %$% 
  InsertRow(.,c("Hamilton", mean(PovertyPercentage)),21) %>% 
  mutate_at(c(2), as.numeric)
  pov$PovertyPercentage <- round(pov$PovertyPercentage,2)

# Retrieve Unemployment Info
une <- data %>% select(CountyName, InLaborForce, Unemployed) %$% 
left_join(setNames(aggregate(InLaborForce, FUN = sum, by = list(CountyName)), 
c("CountyName", "LaborForceSum")), setNames(aggregate(Unemployed, FUN = sum, 
by = list(CountyName)), c("CountyName", "UnemployedSum")), 
by = c("CountyName")) %>% 
mutate(UnemploymentRate = UnemployedSum/LaborForceSum) %$% 
InsertRow(., c("Hamilton",0,0,mean(UnemploymentRate)),21) %>% 
mutate_at(c(2,3,4), as.numeric)

# Create Poverty Boxplot
par(mar = c(1,7,1,1))
Boxplot(pov$PovertyPercentage, id = list(labels = pov$CountyName, cex = 3), ann = FALSE, cex.axis = 2.5, cex = 3)
title(ylab = "Poverty Rate (%)", cex.lab = 3, line = 4.5)

#create Unemployment Boxplot
Boxplot(une$UnemploymentRate, id = list(labels = une$CountyName, cex = 3), ann = FALSE, cex.axis = 2.5, cex = 3)
title(ylab = "Unemployment Rate (%)", cex.lab = 3, line = 4.5)

# MCMC diagnostics
library(coda)
library(ggmcmc)
library(mcmcplots)

setwd("./input_files/coda")
file.List <- list.files(".")
file.List <- file.List[sort.list(file.List)]
# load coda files and generate MCMC objects for each models 
for(i in seq(25, length(file.List), 4)){
  subname <- sub("CODA.*","",file.List[i]) # get the substring before "CODA" 
  chain1 <- read.coda(output.file = paste0(subname, "CODA_chain1.txt"), index.file = paste0(subname, "CODA_index.txt"))
  chain2 <- read.coda(output.file = paste0(subname, "CODA_chain2.txt"), index.file = paste0(subname, "CODA_index.txt"))
  chain3 <- read.coda(output.file = paste0(subname, "CODA_chain3.txt"), index.file = paste0(subname, "CODA_index.txt"))
  var.name <- paste(gsub("[_]","", sub("CODA.*","",file.List[i])), sep = "")
  assign(var.name, as.mcmc.list(list(chain1, chain2, chain3))) # create variable dynamically 
  assign(paste0(var.name, ".gg"), ggs(get(var.name))) #retrieve the variable object as input to ggs function 
}

glo.var <- names(as.list(.GlobalEnv)) %>% .[grep(".gg", .)] %>% 
str_remove_all(., ".gg") # get all variable name that is not with .gg
# to generate all html files with trace, density, autocorrelation plot in one page
for(i in 1:length(glo.var)){
  dir.create(paste0("/Users/jiayan/Downloads/HTML_Diagnostics/",glo.var[i])) #create subfolder to store all html files
  mcmcplot(get(glo.var[i]), style = "plain", greek = T, 
  dir = paste0("/Users/jiayan/Downloads/HTML_Diagnostics/",glo.var[i]),
  filename = glo.var[i], extension = "html", heading = glo.var[i])
}

# Gelman & Rubin diagnostics 
# CAUTION: high computing power, run function without loops if necessary
for(i in 1:length(glo.var)){
  gelman.diag(get(glo.var[i]))
  gelman.plot(get(glo.var[i]))
}

# to generate boxplot for each parameters using the chains value retrieved from coda files 
dc.M4EinterT4E <- M4EinterT4E.gg %>% filter(Parameter %in% c("dc[1]", "dc[2]", "dc[3]", 
"dc[4]", "dc[5]")) %>% select(3,4) %>% mutate_at(c(1), as.character) %>% 
mutate(param = case_when(Parameter == "dc[1]" ~ "2000", 
Parameter == "dc[2]" ~ "2001", Parameter == "dc[3]" ~ "2002", 
Parameter == "dc[4]" ~ "2003", Parameter == "dc[5]" ~ "2004")) %>% 
select(param, value)
par(mar = c(4,5,1,1))
# create M4EinterT4E temporal boxplot 
boxplot(value ~ param, dc.M4EinterT4E, ylab = expression(paste("c"[t])), 
        xlab = "", cex.lab = 2, cex.axis = 2)

dc.M2inter <- M2inter.gg %>% filter(Parameter %in% c("dc[1]", "dc[2]", "dc[3]", 
"dc[4]", "dc[5]")) %>% select(3,4) %>% mutate_at(c(1), as.character) %>% 
mutate(param = case_when(Parameter == "dc[1]" ~ "2000", 
Parameter == "dc[2]" ~ "2001", Parameter == "dc[3]" ~ "2002", 
Parameter == "dc[4]" ~ "2003", Parameter == "dc[5]" ~ "2004")) %>% 
select(param, value)
# create M2inter temporal boxplot 
boxplot(value ~ param, dc.m2inter, ylab = expression(paste("c"[t])), xlab = "", 
        cex.lab = 2, cex.axis = 2)

# to generate boxplot for spatial precision prior, tau
tau.M4EinterT4E <- M4EinterT4E.gg %>% filter(Parameter %in% c("tau[1]", 
"tau[2]", "tau[3]", "tau[4]", "tau[5]")) %>% select(3,4) %>% 
mutate_at(c(1), as.character) %>% 
mutate(param = case_when(Parameter == "tau[1]" ~ "tau^1", 
Parameter == "tau[2]" ~ "tau^2", Parameter == "tau[3]" ~ "tau^3", 
Parameter == "tau[4]" ~ "tau^4", Parameter == "tau[5]" ~ "tau^5")) %>% 
select(param, value)
# create M4EinterT4E spatial precision boxplot
boxplot(value ~ param, tau.M4EinterT4E, ylab = "values", xlab = "", 
        cex.lab = 2, cex.axis = 2, xaxt = "n")
axis(1, at=1:5, labels=c(parse(text=paste("tau^(1)")),
parse(text=paste("tau^(2)")),parse(text=paste("tau^(3)")),
parse(text=paste("tau^(4)")),parse(text=paste("tau^(5)"))), cex.axis = 2)

tau.M2inter <- M2inter.gg %>% filter(Parameter %in% c("tau[1]", "tau[2]", 
"tau[3]", "tau[4]", "tau[5]")) %>% select(3,4) %>% mutate_at(c(1), 
as.character) %>% mutate(param = case_when(Parameter == "tau[1]" ~ "tau^1", 
Parameter == "tau[2]" ~ "tau^2", Parameter == "tau[3]" ~ "tau^3", 
Parameter == "tau[4]" ~ "tau^4", Parameter == "tau[5]" ~ "tau^5")) %>% 
select(param, value)
# create M2inter spatial precision boxplot
boxplot(value ~ param, tau.M2inter, ylab = "values", xlab = "", cex.lab = 2, 
        cex.axis = 2, xaxt = "n")
axis(1, at=1:5, labels=c(parse(text=paste("tau^(1)")),
parse(text=paste("tau^(2)")),parse(text=paste("tau^(3)")),
parse(text=paste("tau^(4)")),parse(text=paste("tau^(5)"))), cex.axis = 2)

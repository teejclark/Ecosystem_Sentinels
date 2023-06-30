
library(tidyverse)
library(googlesheets4)
library(metafor)
library(glmulti)

##################################################################################################################

# import data
df <- read.csv("Ecosystem Sentinels Datasheet.csv")

##################################################################################################################

#### Create a Map ####
library(rnaturalearth)
library(sp)
library(tmap)
library(viridis)
library(patchwork)

# filter data
map_data <- df %>%
  distinct(Paper, .keep_all = T) %>% # keep distinct papers
  filter(!is.na(Longitude)) %>% # filter NAs in long/lat
  filter(!is.na(Latitude)) %>%
  {SpatialPointsDataFrame(cbind(.$Longitude, .$Latitude), .,
                          proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# Briana's CCHWC paper map - using ggplot
# World map with locations of case studies
world <- map_data("world")

a <- df %>%
  distinct(Paper, .keep_all = T) %>% # keep distinct papers
  filter(!is.na(Longitude)) %>% # filter NAs in long/lat
  filter(!is.na(Latitude)) %>%
  ggplot()+
  geom_map(data = world, map = world, aes(x=long, y=lat, map_id = region), fill = "lightgray")+
  geom_point(aes(x=Longitude, y = Latitude, color = Taxonomy, shape = System), size = 5, alpha = 1)+# size = 4, alpha = 0.85
  theme_void()+
  scale_color_viridis(option = "turbo", direction = 1, discrete = TRUE, name = "Taxonomy")+
  guides(shape = guide_legend(title = "System"))+
  theme(legend.text = element_text(size = 16), legend.title = element_text(size=24))

# Histogram of Sentinel Measurement
b <- df %>%
  filter(`Sampling Method` != "NA") %>%
  ggplot(aes(x = `Sampling Method`, fill = System))+
  geom_histogram(stat = "count", color = "black", size = 1)+
  scale_fill_manual(values = c("#ffd166","#ef476f","#06d6a0"))+
  theme_classic()+
  ylab("# samples")+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=24),
    axis.text = element_text(size=18, color="black"),
    legend.text = element_text(size=16,color="black"),
    legend.position = "none",
    axis.text.x = element_text(angle = 50, hjust=1))

# Histogram of Environmental Category
c <- df %>%
  filter(`Connection? A or B?` == "B") %>%
  ggplot(aes(x = `Environmental Category`, fill = System))+
  geom_histogram(stat = "count", color = "black", size = 1)+
  scale_fill_manual(values = c("#ffd166","#ef476f","#06d6a0"))+
  theme_classic()+
  ylab("# samples")+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=24),
    axis.text = element_text(size=18, color="black"),
    legend.text = element_text(size=16,color="black"),
    legend.position = "none",
    axis.text.x = element_text(angle = 50, hjust=1))

# Histogram of Ecosystem Category
d <- df %>%
  filter(`Connection? A or B?` == "A") %>%
  ggplot(aes(x = `Ecosystem Category`, fill = System))+
  geom_histogram(stat = "count", color = "black", size = 1)+
  scale_fill_manual(values = c("#ffd166","#ef476f","#06d6a0"))+
  theme_classic()+
  ylab("# samples")+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=24),
    axis.text = element_text(size=18, color="black"),
    legend.text = element_text(size=16,color="black"),
    legend.title = element_text(size = 24, color = "black"),
    axis.text.x = element_text(angle = 50, hjust=1))

# Histogram of specific taxonomy
e <- df %>%
  distinct(Paper, .keep_all = T) %>% # keep distinct papers
  filter(`Revised Taxnomy` != "NA") %>%
  ggplot(aes(x = `Revised Taxnomy`, fill = System)) +
  geom_histogram(stat = "count", color = "black", size = 0.5)+
  theme_classic()+
  coord_flip()


##################################################################################################################

#### Basic Meta-analysis Calculations ####

# first calculate the transformed data
# removes NAs
blah <- df %>% filter(Value != "NA")

# NOTE: need to calculate on absolute value? right?
dat <- escalc(measure = "COR", ri = abs(Value), ni = `Sample Size`,
              data = blah, vtype = "AV")

# calculate a sample meta-analysis, investigate impact of weights
dat$ni <- dat$Sample.Size

# assign unique ID to each study (should be 114 unique studies)
dat2 <- dat %>%
  group_by(Paper) %>%
  mutate(ID = cur_group_id())

# initial meta-analysis without variables
res5 <- rma.mv(yi, vi, W = log(ni), data = dat2, random = list(~1|ID))

# let's do backwards stepwise model selection

res_full <- rma.mv(yi, vi, W = log(ni), 
                  mods = ~Taxonomy + System + Trophic.Position + Ecto.vs..Endotherm. + Connection..A.or.B. +
                  + Sampling.Method + as.numeric(Body.Size..kg.) + Scale.of.Change + 
                    Environmental.Category + Ecosystem.Category,
               data = dat2, random = list(~1|ID))

res_sel <- rma.mv(yi, vi, W = log(ni), 
                   mods = ~Connection..A.or.B. + Trophic.Position + Ecto.vs..Endotherm. + 
                    Scale.of.Change + Environmental.Category,
                   data = dat2, random = list(~1|ID)); AIC(res_sel) # FINAL SELECTION!!!

# final variables: connection, trophic position, ecto vs. endotherm, scale of change, environmental category

# testing A and B connections independently
dat2_A <- dat2 %>% filter(Connection..A.or.B. == "A")

res_sel <- rma.mv(yi, vi, W = log(ni), 
                   mods = ~Taxonomy + System + Trophic.Position + Ecto.vs..Endotherm. + 
                     + Sampling.Method + as.numeric(Body.Size..kg.) + Scale.of.Change + 
                     Environmental.Category + Ecosystem.Category,
                   data = dat2_A, random = list(~1|ID)); AIC(res_sel)

res_final <- rma.mv(yi, vi, W = log(ni), 
                  mods = ~Trophic.Position + as.numeric(Body.Size..kg.) + Ecosystem.Category,
                  data = dat2_A, random = list(~1|ID))

dat2_B <- dat2 %>% filter(Connection..A.or.B. == "B")

res_sel <- rma.mv(yi, vi, W = log(ni), 
                  mods = ~Taxonomy + System + Trophic.Position + Ecto.vs..Endotherm. + 
                    + Sampling.Method + as.numeric(Body.Size..kg.) + Scale.of.Change + 
                    Environmental.Category + Ecosystem.Category,
                  data = dat2_B, random = list(~1|ID)); AIC(res_sel)

res_final <- rma.mv(yi, vi, W = log(ni), 
                  mods = ~Trophic.Position + Environmental.Category,
                  data = dat2_B, random = list(~1|ID))

##################################################################################################################

#### OTHER GRAPHS ####

# Sampling method - effect size
df %>%
  filter(`Sampling Method` != "NA") %>%
  ggplot(aes(x = `Sampling Method`, y = abs(Value)))+
  geom_boxplot(fill = "#E5E255")+
  theme_classic()+
  ylab("Effect Size")+
  theme(text = element_text(size = 12),
        axis.text = element_text(color = "black"))

# Sampling method - vote counting
blah <- df %>%
  group_by(`Sampling Method`,`Vote Counting`) %>%
  summarize(total = n()) %>%
  mutate(freq = round(total/sum(total),3))

blah %>%
  filter(`Sampling Method` != "NA") %>%
  ggplot(aes(x = `Sampling Method`, y=freq, fill = as.factor(`Vote Counting`)))+
  geom_col()

# Environmental Category - effect size
df %>%
  filter(`Connection? A or B?` == "B") %>%
  ggplot(aes(x = `Environmental Category`, y = abs(Value)))+
  geom_boxplot(fill = "#EC6D2F")+
  theme_classic()+
  ylab("Effect Size")+
  theme(text = element_text(size = 12),
        axis.text = element_text(color = "black"))
  
# Ecosystem Category - effect size
df %>%
  filter(`Connection? A or B?` == "A") %>%
  ggplot(aes(x = `Ecosystem Category`, y = abs(Value)))+
  geom_boxplot(fill = "#8E2929")+
  theme_classic()+
  ylab("Effect Size")+
  theme(text = element_text(size = 12),
        axis.text = element_text(color = "black"))

#######################################################################################################

#### Final Summary Figures

# combine trophic position and endo/ecto
df$`Trophic Position`[df$`Trophic Position` == "Primary Consumer"] <- "1° Consumer"
df$`Trophic Position`[df$`Trophic Position` == "Secondary Consumer"] <- "2° Consumer"

a <- df %>%
  ggplot(aes(x = `Trophic Position`, y = abs(Value)))+
  geom_boxplot(fill = "#756FB4")+
  theme_classic()+
  ylab("Effect Size")+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=24),
    axis.text = element_text(size=18, color="black"))

b <- df %>%
  filter(`Connection? A or B?` == "B") %>%
  ggplot(aes(x = `Environmental Category`, y = abs(Value)))+
  geom_boxplot(fill = "#D96228")+
  theme_classic()+
  ylab("Effect Size")+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=24),
    axis.text = element_text(size=18, color="black"),
    axis.text.x = element_text(angle = 30, hjust=1))

# supplementary figure
c <- df %>%
  filter(`Connection? A or B?` == "A") %>%
  ggplot(aes(x = `Ecosystem Category`, y = abs(Value)))+
  geom_boxplot(fill = "#1F9F78")+
  theme_classic()+
  ylab("Effect Size")+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=24),
    axis.text = element_text(size=18, color="black"))

c <- df %>%
  ggplot(aes(x = `Ecto vs. Endotherm?`, y = abs(Value)))+
  geom_boxplot(fatten = NULL, fill = "#DBFCFF")+
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid")+
  theme_classic()+
  xlab("Thermoregulation")+
  ylab("Effect Size")+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=24),
    axis.text = element_text(size=18, color="black"))

d <- df %>%
  filter(`Scale of Change` != "Hours") %>%
  filter(`Scale of Change` != "Months") %>%
  ggplot(aes(x = `Scale of Change`, y = abs(Value)))+
  geom_boxplot(fill = "#FFF275")+
  theme_classic()+
  xlab("Sampling Scale")+
  ylab("Effect Size")+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=24),
    axis.text = element_text(size=18, color="black"))

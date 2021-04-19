##############################
## ---- Working directory ----
##############################

setwd(dir = '/Users/s2018147/Documents/Working directory/Mouse')
source(file = '/Users/s2018147/Documents/Working directory/gene_drive_mouse.R')

#############################
## ---- Parameter inputs ----
#############################

input <- list()

input$strategy <- 5 ## 1 = Heterozygotic XX sterility
                    ## 2 = Heterozygotic XX sex reversal
                    ## 3 = Homozygotic embryonic non-viability
                    ## 4 = Homozygotic XX sterility
                    ## 5 = Heterozygotic XX sterility
input$nIter <- 100 # number of iterations
input$nYear <- 20 # maximum number of years
input$K <- 50000 # carrying capacity of offspring
input$nG <-  100 # initial number carrying the gene drive
input$supplement <- seq(0,1,by=0.05) # amount of animals supplemented per supplementation
input$suppInterval <- seq(0.1,1,by=0.1) # interval at which supplementation happens in years
input$nGuides <- 3 # number of guide RNAs

input$Pc <- 0.95 # cutting probability
input$Pn <- 0.02 # NHEJ probability
input$m <- 6 # litter size
input$genTime <- 5.2 # generation time (in weeks)
input$rmax <- 7.76 # maximum annual population growth rate
input$Pr <- 0 # probability of inherent resistance (affects number of susceptible sites at initiation)
input$pRnonFunc <- 2/3 # probability of single-site deletion causing loss-of-function (if position=='exon')
input$position <- 'exon' # c('intron', 'exon') # position of gene drive, one of 'intron' or 'exon'
input$transType <- 'sim' # c('seq', 'sim') # transition matrix type, either 'seq' (sequential) or 'sim' (simultaneous)
input$compensation <- FALSE ## reproductive compensation, for strategy 3 only
input$implantSmax <- 1 # maximum embryonic survival rate (when compensation==TRUE)
input$implantSmin <- 1 # minimum embryonic survival rate (when compensation==TRUE)
input$implantTheta <- 1 # theta parameter for theta-logistic compensation (when compensation==TRUE)

inputs <- expand.grid(input)
View(inputs)
rownames(inputs) <- 1:nrow(inputs)

#######################################
## ---- Set up parallel processing ----
#######################################

nproc <- 7
cl.tmp = makeCluster(rep('localhost',nproc), type = 'SOCK')
registerDoSNOW(cl.tmp)
getDoParWorkers()

############################
## ---- Run simulations ----
############################

## run just one
# res <- gd.comp(input = inputs[2, ])

## run in parallel
res <- foreach(rowNum = 1:nrow(inputs), .packages = 'reshape2', .verbose = TRUE) %dopar% { gd.comp(input = inputs[rowNum, ]) }

# make a nice form for further processing
res <- do.call(rbind, res)
var <- c('nInd', 'nIndF', 'nIndM',
         'nIndWW', 'nIndWG', 'nIndWN', 'nIndGG', 'nIndGN', 'nIndNN',
         'pIndWW', 'pIndWG', 'pIndWN', 'pIndGG', 'pIndGN', 'pIndNN',
         'nIndWWF', 'nIndWGF', 'nIndWNF', 'nIndGGF', 'nIndGNF', 'nIndNNF',
         'pIndWWF', 'pIndWGF', 'pIndWNF', 'pIndGGF', 'pIndGNF', 'pIndNNF',
         'nIndWWM', 'nIndWGM', 'nIndWNM', 'nIndGGM', 'nIndGNM', 'nIndNNM',
         'pIndWWM', 'pIndWGM', 'pIndWNM', 'pIndGGM', 'pIndGNM', 'pIndNNM',
         'nAllW', 'nAllG', 'nAllN',
         'pAllW', 'pAllG', 'pAllN',
         'nAllWF', 'nAllGF', 'nAllNF',
         'pAllWF', 'pAllGF', 'pAllNF',
         'nAllWM', 'nAllGM', 'nAllNM',
         'pAllWM', 'pAllGM', 'pAllNM')
res[, var] = lapply(X = res[, var], FUN = function(z) { z[is.na(z)] = 0; z })
res = as_tibble(res)

############################
## ---- Process results ----
############################

## collate and calculate population parameters over time
res$strategy <- factor(res$strategy)
res$supplement <- factor(res$supplement)
res$suppInterval <- factor(res$suppInterval)

res = res %>%
  filter(year > 0) %>%
  mutate(extinct = as.numeric(nInd == 0) * 100,
         pInd = nInd / K * 100)

resSum = res %>%
  group_by(strategy, supplement, suppInterval, year) %>%
  summarise(extinctMean = mean(extinct, na.rm = TRUE),
            extinctQL   = quantile(extinct, na.rm = TRUE, probs = 0.025),
            extinctQU   = quantile(extinct, na.rm = TRUE, probs = 0.975),

            nIndMean = mean(nInd, na.rm = TRUE),
            nIndQL   = quantile(nInd, na.rm = TRUE, probs = 0.025),
            nIndQU   = quantile(nInd, na.rm = TRUE, probs = 0.975),

            pIndMean = mean(pInd, na.rm = TRUE),
            pIndQL   = quantile(pInd, na.rm = TRUE, probs = 0.025),
            pIndQU   = quantile(pInd, na.rm = TRUE, probs = 0.975),

            pAllWMean = mean(pAllW, na.rm = TRUE),
            pAllWQL   = quantile(pAllW, na.rm = TRUE, probs = 0.025),
            pAllWQU   = quantile(pAllW, na.rm = TRUE, probs = 0.975),

            pAllGMean = mean(pAllG, na.rm = TRUE),
            pAllGQL   = quantile(pAllG, na.rm = TRUE, probs = 0.025),
            pAllGQU   = quantile(pAllG, na.rm = TRUE, probs = 0.975),

            pAllNMean = mean(pAllN, na.rm = TRUE),
            pAllNQL   = quantile(pAllN, na.rm = TRUE, probs = 0.025),
            pAllNQU   = quantile(pAllN, na.rm = TRUE, probs = 0.975),

            pAllWFMean = mean(pAllWF, na.rm = TRUE),
            pAllWFQL   = quantile(pAllWF, na.rm = TRUE, probs = 0.025),
            pAllWFQU   = quantile(pAllWF, na.rm = TRUE, probs = 0.975),

            pAllGFMean = mean(pAllGF, na.rm = TRUE),
            pAllGFQL   = quantile(pAllGF, na.rm = TRUE, probs = 0.025),
            pAllGFQU   = quantile(pAllGF, na.rm = TRUE, probs = 0.975),

            pAllNFMean = mean(pAllNF, na.rm = TRUE),
            pAllNFQL   = quantile(pAllNF, na.rm = TRUE, probs = 0.025),
            pAllNFQU   = quantile(pAllNF, na.rm = TRUE, probs = 0.975),

            pAllWMMean = mean(pAllWM, na.rm = TRUE),
            pAllWMQL   = quantile(pAllWM, na.rm = TRUE, probs = 0.025),
            pAllWMQU   = quantile(pAllWM, na.rm = TRUE, probs = 0.975),

            pAllGMMean = mean(pAllGM, na.rm = TRUE),
            pAllGMQL   = quantile(pAllGM, na.rm = TRUE, probs = 0.025),
            pAllGMQU   = quantile(pAllGM, na.rm = TRUE, probs = 0.975),

            pAllNMMean = mean(pAllNM, na.rm = TRUE),
            pAllNMQL   = quantile(pAllNM, na.rm = TRUE, probs = 0.025),
            pAllNMQU   = quantile(pAllNM, na.rm = TRUE, probs = 0.975))

# save(res, resSum, file = 'mouse_strategy5_heatmaps.RData')
# load(file = 'mouse_strategy5_heatmaps.RData')

##################
## ---- Plots ----
##################

# ---- Plot themes ----

PaperTheme <- theme_bw(base_size = 11, base_family = "serif") + theme(legend.position = "right", strip.background = element_blank())
PaperThemeNoLegend <- PaperTheme + theme(legend.position = "none")
PreseTheme <- theme_bw(base_size = 18, base_family = "sans")  + theme(legend.position = "right", strip.background = element_blank())
PreseThemeNoLegend <- PreseTheme + theme(legend.position = "none")

PaperSize <- 10
PreseSize <- 16

# ---- Plotting ----

heatMapData5 <- resSum[resSum$year==5,]
heatMapData10 <- resSum[resSum$year==10,]
heatMapData20 <- resSum[resSum$year==20,]
heatMapData5$extinctRate <- heatMapData5$extinctMean/100
heatMapData10$extinctRate <- heatMapData10$extinctMean/100
heatMapData20$extinctRate <- heatMapData20$extinctMean/100

# Suppression rate after 5 years
p = heatMapData5 %>%
  ggplot() +
  geom_tile(mapping = aes(x = supplement, y = suppInterval, fill = extinctRate)) +
  scale_fill_viridis(option="viridis", name = "Suppresion rate") +
  theme_bw() +
  theme(legend.position="right", strip.background = element_blank()) +
  scale_x_discrete(breaks = seq(0,2000,200)) +
  xlab(label = "Supplementation amount") +
  ylab(label = "Supplementation interval (years)")
p

ggsave(plot = p + PaperTheme, filename = "mouse_fig11_suppression_5yrs.png",    height = PaperSize, width = PaperSize * 1.5, unit = "cm")
ggsave(plot = p + PreseTheme, filename = "mouse_fig11_suppression_5yrs_PPT.png", height = PreseSize, width = PreseSize * 1.5, unit = "cm")

# Suppression rate after 10 years
p = heatMapData10 %>%
  ggplot() +
  geom_tile(mapping = aes(x = supplement, y = suppInterval, fill = extinctRate)) +
  scale_fill_viridis(option="viridis", name = "Suppresion rate") +
  theme_bw() +
  theme(legend.position="right", strip.background = element_blank()) +
  scale_x_discrete(breaks = seq(0,2000,200)) +
  xlab(label = "Supplementation amount") +
  ylab(label = "Supplementation interval (years)")
p

ggsave(plot = p + PaperTheme, filename = "mouse_fig12_suppression_10yrs.png",    height = PaperSize, width = PaperSize * 1.5, unit = "cm")
ggsave(plot = p + PreseTheme, filename = "mouse_fig12_suppression_10yrs_PPT.png", height = PreseSize, width = PreseSize * 1.5, unit = "cm")

# Suppression rate after 20 years
p = heatMapData20 %>%
  ggplot() +
  geom_tile(mapping = aes(x = supplement, y = suppInterval, fill = extinctRate)) +
  scale_fill_viridis(option="viridis", name = "Suppresion rate") +
  theme_bw() +
  theme(legend.position="right", strip.background = element_blank()) +
  scale_x_discrete(breaks = seq(0,2000,200)) +
  xlab(label = "Supplementation amount") +
  ylab(label = "Supplementation interval (years)")
p

ggsave(plot = p + PaperTheme, filename = "mouse_fig13_suppression_20yrs.png",    height = PaperSize, width = PaperSize * 1.5, unit = "cm")
ggsave(plot = p + PreseTheme, filename = "mouse_fig13_suppression_20yrs_PPT.png", height = PreseSize, width = PreseSize * 1.5, unit = "cm")

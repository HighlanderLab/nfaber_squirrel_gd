##############################
## ---- Working directory ----
##############################

setwd(dir = '/Users/s2018147/Documents/Working directory/Squirrel')
source(file = '/Users/s2018147/Documents/Working directory/gene_drive_squirrel.R')

#############################
## ---- Parameter inputs ----
#############################

input <- list()
## 1 = Heterozygotic XX sterility
## 2 = Heterozygotic XX sex reversal
## 3 = Homozygotic embryonic non-viability
## 4 = Homozygotic XX sterility
## 5 = Roslin version of 4
## 6 = X-shredder
## 7 = Homozygotic XX sterility with TA
## 8 = X-shredder with TA

input$strategy <- 7

input$nIter <- 100 # number of iterations
input$nYear <-  50 # maximum number of years

input$K <- 3000 # carrying capacity of offspring
input$nG <-  100 # initial number carrying the gene drive
input$supplement <- 0 # c(0,25,50,100,200,500) # amount of animals supplemented per supplementation
input$suppInterval <- 0 # c(1,2,3,4,5) # interval at which supplementation happens in years

input$nGuides <- 1 # number of guide RNA target sites
input$Pc <- 0.95 # cutting probability
input$Pn <- c(0.02,0.1,0.5) # NHEJ probability
input$daisyfield <- -1 # c(10,20,40,60,100,500) # number of daisies in daisy field, -1 if no daisy drive. 
                                      # Must be an even number as daisy sites will be homozygous
input$X.shr.eff <- 0 # c(0.5,0.6,0.7,0.8,0.9,1) # percent of male offspring when drive present in fathers
input$nGuides.TA <- c(1,2,4)
input$overexpressionOK <- c(TRUE,FALSE)
input$lethal <- c(TRUE,FALSE)
input$m <- 2.87 # litter size
input$genTime <- 26 # generation time (in weeks)
input$rmax <- log(3.2) # maximum annual population growth rate
input$Pr <- 0 # probability of inherent resistance (affects number of susceptible sites at initiation)
input$pRnonFunc <- 2/3 # probability of single-site deletion causing loss-of-function (if position=='exon')
input$pDoubleNHEJ <- 0.2 # probability that a block is not removed when two gRNAs cut on the outside of it
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

## run one
# res <- gd.comp(input = inputs[16, ])

## run in parallel
res <- foreach(rowNum = 1:nrow(inputs), .packages = 'reshape2', .verbose = TRUE) %dopar% { gd.comp(input = inputs[rowNum, ]) }

# make a nice form for further processing
res <- do.call(rbind, res)
# View(head(res, n=100))
var <- c('nInd', 'nIndF', 'nIndM',
         'nIndWW', 'nIndWG', 'nIndWN', 'nIndGG', 'nIndGN', 'nIndNN',
         'pIndWW', 'pIndWG', 'pIndWN', 'pIndGG', 'pIndGN', 'pIndNN',
         'nIndWWF', 'nIndWGF', 'nIndWNF', 'nIndGGF', 'nIndGNF', 'nIndNNF',
         'pIndWWF', 'pIndWGF', 'pIndWNF', 'pIndGGF', 'pIndGNF', 'pIndNNF',
         'nIndWWM', 'nIndWGM', 'nIndWNM', 'nIndGGM', 'nIndGNM', 'nIndNNM',
         'pIndWWM', 'pIndWGM', 'pIndWNM', 'pIndGGM', 'pIndGNM', 'pIndNNM',
         'nAllW', 'nAllG', 'nAllN',
         'pAllW', 'pAllG', 'pAllN',
         'pAllW.TA','pAllG.TA','pAllN.TA',
         'nAllWF', 'nAllGF', 'nAllNF',
         'pAllWF', 'pAllGF', 'pAllNF',
         'nAllWM', 'nAllGM', 'nAllNM',
         'pAllWM', 'pAllGM', 'pAllNM')
res[, var] = lapply(X = res[, var], FUN = function(z) { z[is.na(z)] = 0; z })
res = as_tibble(res)
# save(res, file = 'mouse_plots_data_raw_bottom.RData')
# View(head(res, n=100))

############################
## ---- Process results ----
############################

## collate and calculate population parameters over time
res$strategy <- factor(res$strategy)
res$Pn <- factor(res$Pn)
res$nGuides <- factor(res$nGuides)
res$supplement <- factor(res$supplement)
res$X.shr.eff <- factor(res$X.shr.eff)
res$daisyfield <- factor(res$daisyfield)
res$suppInterval <- factor(res$suppInterval)
res$nGuides.TA <- factor(res$nGuides.TA)
res$lethal <- factor(res$lethal)
res$overexpressionOK <- factor(res$overexpressionOK)

res = res %>%
  mutate(extinct = as.numeric(nInd == 0) * 100,
         pInd = nInd / K * 100)

resSum = res %>%
  group_by(strategy,nGuides,supplement,suppInterval,daisyfield,Pn,X.shr.eff,nGuides.TA,lethal,overexpressionOK,year) %>%
  summarise(extinctMean = mean(extinct, na.rm = TRUE),
            extinctQL   = quantile(extinct, na.rm = TRUE, probs = 0.025),
            extinctQU   = quantile(extinct, na.rm = TRUE, probs = 0.975),

            nIndMean = mean(nInd, na.rm = TRUE),
            nIndQL   = quantile(nInd, na.rm = TRUE, probs = 0.025),
            nIndQU   = quantile(nInd, na.rm = TRUE, probs = 0.975),
            
            nMaleMean = mean(nIndM, na.rm = TRUE),
            nMaleQL   = quantile(nIndM, na.rm = TRUE, probs = 0.025),
            nMaleQU   = quantile(nIndM, na.rm = TRUE, probs = 0.975),
            
            nFemaleMean = mean(nIndF, na.rm = TRUE),
            nFemaleQL   = quantile(nIndF, na.rm = TRUE, probs = 0.025),
            nFemaleQU   = quantile(nIndF, na.rm = TRUE, probs = 0.975),
            
            male.femaleRatio = mean(ifelse(nIndF+nIndM==0,0,nIndF/(nIndF+nIndM)), na.rm = TRUE),
            male.femaleRatioQL   = quantile(ifelse(nIndF+nIndM==0,0,nIndF/(nIndF+nIndM)), na.rm = TRUE, probs = 0.025),
            male.femaleRatioQU   = quantile(ifelse(nIndF+nIndM==0,0,nIndF/(nIndF+nIndM)), na.rm = TRUE, probs = 0.975),

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
            
            pAllW.TAMean = mean(pAllW.TA, na.rm = TRUE),
            pAllW.TAQL = quantile(pAllW.TA, na.rm = TRUE, probs = 0.025),
            pAllW.TAQU = quantile(pAllW.TA, na.rm = TRUE, probs = 0.975),
            
            pAllG.TAMean = mean(pAllG.TA, na.rm = TRUE),
            pAllG.TAQL = quantile(pAllG.TA, na.rm = TRUE, probs = 0.025),
            pAllG.TAQU = quantile(pAllG.TA, na.rm = TRUE, probs = 0.975),
            
            pAllN.TAMean = mean(pAllN.TA, na.rm = TRUE),
            pAllN.TAQL = quantile(pAllN.TA, na.rm = TRUE, probs = 0.025),
            pAllN.TAQU = quantile(pAllN.TA, na.rm = TRUE, probs = 0.975),

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

# View(resSum)

# save(res, resSum, file = 'Fig2_supp_whichgenes.RData')
# load(file = 'Fig2_supp_whichgenes.RData')

##################
## ---- Plots ----
##################

# ---- Plot themes ----

PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + theme(legend.position = "bottom", strip.background = element_blank(),panel.grid = element_blank(), legend.title.align = 0.5)
PreseTheme <- theme_bw(base_size = 18, base_family = "sans")  + theme(legend.position = "bottom", strip.background = element_blank(),panel.grid = element_blank(), legend.title.align = 0.5)

PaperSize <- 10
PreseSize <- 16

# ---- Plotting ----

copydf <- resSum
copydf <- rename(copydf, gRNAs = nGuides.TA)
copydf$lethal <- ifelse(as.logical(levels(copydf$lethal)[copydf$lethal]), "Embryonic lethal", "Both-sex infertility")
copydf$overexpressionOK <- ifelse(as.logical(levels(copydf$overexpressionOK)[copydf$overexpressionOK]), "Allowed", "Not allowed")
copydf$lethal <- factor(copydf$lethal, levels=c("Embryonic lethal", "Both-sex infertility"))
copydf$overexpressionOK <- factor(copydf$overexpressionOK, levels=c("Not allowed", "Allowed"))

# nGuides Census
p = copydf %>%
  ggplot() +
  geom_line(mapping = aes(x = year, y = nIndMean, colour = lethal:overexpressionOK)) +
  geom_ribbon(mapping = aes(x = year, ymin = nIndQL, ymax = nIndQU, fill = lethal:overexpressionOK), alpha = .1) +
  scale_colour_viridis(option="viridis", discrete=TRUE, name = "Target gene : Overexpression",guide=guide_legend(ncol=1, title.position = "top")) +
  scale_fill_viridis(option="viridis", discrete=TRUE, name = "Target gene : Overexpression",guide=guide_legend(ncol=1, title.position = "top")) +
  facet_grid(gRNAs ~ Pn, labeller = label_both) +
  theme_bw() +
  theme(legend.position="right", strip.background = element_blank()) +
  ylim(c(0, NA)) +
  xlab(label = "Year") +
  ylab(label = "Population size")
p

ggsave(plot = p + PaperTheme, filename = "Fig2_supp_whichgenes.png",    height = PaperSize * 1.5, width = PaperSize * 1.5, unit = "cm")
ggsave(plot = p + PreseTheme, filename = "Fig2_supp_whichgenes_PPT.png", height = PreseSize * 1.5, width = PreseSize * 1.5, unit = "cm")

##############################
## ---- Working directory ----
##############################

setwd(dir = '/Users/s2018147/Documents/Working directory/Mouse')
source(file = '/Users/s2018147/Documents/Working directory/gene_drive_mouse.R')

#############################
## ---- Parameter inputs ----
#############################

input <- list()

input$strategy <- 4 ## 1 = Heterozygotic XX sterility
                    ## 2 = Heterozygotic XX sex reversal
                    ## 3 = Homozygotic embryonic non-viability
                    ## 4 = Homozygotic XX sterility
                    ## 5 = Heterozygotic XX sterility
input$nIter <- 3 # number of iterations
input$nYear <- 10 # maximum number of years
input$K <- 50000 # carrying capacity of offspring
input$nG <- 100 # initial number carrying the gene drive
input$supplement <- 0 # amount of animals supplemented per supplementation
input$suppInterval <- 0 # interval at which supplementation happens in years
input$nGuides <- c(1,2,3,4,5) # number of guide RNAs

input$Pc <- 0.95 # cutting probability
input$Pn <- c(0,0.02,0.1) # NHEJ probability
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
# res <- gd.comp(input = inputs[1, ])

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
res$Pn <- factor(res$Pn)
res$nGuides <- factor(res$nGuides)

res = res %>%
  filter(year > 0) %>%
  mutate(extinct = as.numeric(nInd == 0) * 100,
         pInd = nInd / K * 100)

resSum = res %>%
  group_by(strategy, nGuides, Pn, year) %>%
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

# save(res, resSum, file = 'mouse_strategy4.RData')
# load(file = 'mouse_strategy4.RData')

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

# nGuides Census
p = resSum %>%
  ggplot() +
  geom_line(mapping = aes(x = year, y = nIndMean, colour = nGuides)) +
  geom_ribbon(mapping = aes(x = year, ymin = nIndQL, ymax = nIndQU, fill = nGuides), alpha = .2) +
  scale_colour_viridis(option="viridis", discrete=TRUE, name = "gRNAs") +
  scale_fill_viridis(option="viridis", discrete=TRUE, name = "gRNAs") +
  facet_grid(. ~ Pn, labeller = label_both) +
  theme_bw() +
  theme(legend.position="right", strip.background = element_blank()) +
  ylim(c(0, NA)) +
  xlab(label = "Year") +
  ylab(label = "Population size")
p

ggsave(plot = p + PaperTheme, filename = "mouse_fig1_populationsize.png",    height = PaperSize, width = PaperSize * 3, unit = "cm")
ggsave(plot = p + PreseTheme, filename = "mouse_fig1_populationsize_PPT.png", height = PreseSize, width = PreseSize * 3, unit = "cm")

# Extinction
p = resSum %>%
  ggplot() +
  geom_line(mapping = aes(x = year, y = extinctMean, colour = nGuides)) +
  scale_colour_viridis(option="viridis", discrete=TRUE, name = "gRNAs") +
  scale_fill_viridis(option="viridis", discrete=TRUE, name = "gRNAs") +
  facet_grid( ~ Pn, labeller = label_both) +
  theme_bw() +
  theme(legend.position="right", strip.background = element_blank()) +
  xlab(label = "Year") +
  ylab(label = "% Eradicated populations")
p

ggsave(plot = p + PaperTheme, filename = "mouse_fig2_suppression.png",    height = PaperSize, width = PaperSize * 3, unit = "cm")
ggsave(plot = p + PreseTheme, filename = "mouse_fig2_suppression_PPT.png", height = PreseSize, width = PreseSize * 3, unit = "cm")

# Allele frequencies
p = resSum %>%
  ggplot() +
  geom_line(mapping = aes(x = year, y = pAllWMean, colour = nGuides)) +
  geom_ribbon(mapping = aes(x = year, ymin = pAllWQL, ymax = pAllWQU, fill = nGuides), alpha = .2) +
  scale_colour_viridis(option="viridis", discrete=TRUE, name = "gRNAs") +
  scale_fill_viridis(option="viridis", discrete=TRUE, name = "gRNAs") +
  facet_grid( ~ Pn, labeller = label_both) +
  theme_bw() +
  theme(legend.position="right", strip.background = element_blank()) +
  xlab(label = "Year") +
  ylab(label = "Wildtype allele frequency") +
  ylim(c(0,1))
p

ggsave(plot = p + PaperTheme, filename = "mouse_fig3_allelewildtype.png",    height = PaperSize, width = PaperSize * 3, unit = "cm")
ggsave(plot = p + PreseTheme, filename = "mouse_fig3_allelewildtype_PPT.png", height = PreseSize, width = PreseSize * 3, unit = "cm")

p = resSum %>%
  ggplot() +
  geom_line(mapping = aes(x = year, y = pAllGMean, colour = nGuides)) +
  geom_ribbon(mapping = aes(x = year, ymin = pAllGQL, ymax = pAllGQU, fill = nGuides), alpha = .2) +
  scale_colour_viridis(option="viridis", discrete=TRUE, name = "gRNAs") +
  scale_fill_viridis(option="viridis", discrete=TRUE, name = "gRNAs") +
  facet_grid( ~ Pn, labeller = label_both) +
  theme_bw() +
  theme(legend.position="right", strip.background = element_blank()) +
  xlab(label = "Year") +
  ylab(label = "Gene drive allele frequency") +
  ylim(c(0,1))
p

ggsave(plot = p + PaperTheme, filename = "mouse_fig4_allelegenedrive.png",    height = PaperSize, width = PaperSize * 3, unit = "cm")
ggsave(plot = p + PreseTheme, filename = "mouse_fig4_allelegenedrive_PPT.png", height = PreseSize, width = PreseSize * 3, unit = "cm")

p = resSum %>%
  ggplot() +
  geom_line(mapping = aes(x = year, y = pAllNMean, colour = nGuides)) +
  geom_ribbon(mapping = aes(x = year, ymin = pAllNQL, ymax = pAllNQU, fill = nGuides), alpha = .2) +
  scale_colour_viridis(option="viridis", discrete=TRUE, name = "gRNAs") +
  scale_fill_viridis(option="viridis", discrete=TRUE, name = "gRNAs") +
  facet_grid( ~ Pn, labeller = label_both) +
  theme_bw() +
  theme(legend.position="right", strip.background = element_blank()) +
  xlab(label = "Year") +
  ylab(label = "Resistant allele frequency") +
  ylim(c(0,1))
p

ggsave(plot = p + PaperTheme, filename = "mouse_fig5_alleleresistant.png",    height = PaperSize, width = PaperSize * 3, unit = "cm")
ggsave(plot = p + PreseTheme, filename = "mouse_fig5_alleleresistant_PPT.png", height = PreseSize, width = PreseSize * 3, unit = "cm")

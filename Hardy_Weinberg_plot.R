library(tidyverse)


########
## Hardy-Weinberg Code by: Mark Ravinet & Glenn-Peter Sætre Oslo, October 2018
## https://evolutionarygenetics.github.io/Chapter3.html

########

p <- 0.8
q <- 1 - p 
# check p and q are equal to 1
(q + p) == 1
# calculate the expected genotype frequencies (_e denotes expected)
A1A1_e <- p^2
A1A2_e <- 2 * (p * q)
A2A2_e <- q^2
# show the allele frequencies in the console
c(A1A1_e, A1A2_e, A2A2_e)
# since these are genotype frequencies, they should also sum to 1 - you can check this like so
sum(c(A1A1_e, A1A2_e, A2A2_e))

# generate a range for p
p <- seq(0, 1, 0.01)
# and also for q
q <- 1 - p

# generate the expected genotype frequencies
A1A1 <- p^2
A1A2 <- 2 * (p * q)
A2A2 <- q^2

# arrange allele frequencies into a tibble/data.frame
geno_freq <- as.tibble(cbind(p, q, A1A1, A1A2, A2A2))

# Use gather to reshape the data.frame for straightforward plotting
geno_freq <- gather(geno_freq, key = "genotype", value = "freq", -p, -q)

# plot the expected genotype frequencies
a <- ggplot(geno_freq, aes(p, freq, colour = genotype)) + geom_line()
a <-a + ylab("Genotype frequency") + xlab("Gene frequency of A1\n(p frequency)")
a + theme_light() + theme(legend.position = "bottom")









########
# All following code written by: Joey White
# for CSCI 2897 final project
########

#Inbreeding effects
gens <- c(0:10)


H <- (1/(2^gens))


curve((1/(2^x)), from=0, to=10, main = "Decrease in Heterozygosity per \nGeneration of Inbreeding", 
      xlab="Generation (Fx)", 
      ylab="Heterozygote Porbability",
      col = "FireBrick",
      lwd = 2)
points (H~gens)




### H and P+Q proportions per generation 

P_Q = 1-H

curve((1/(2^x)), from=0, to=10, main = "Genotype Proportions per \nGeneration of Inbreeding", 
      xlab="Generation (Fx)", 
      ylab="Proportion",
      col = "FireBrick",
      lwd = 2)
curve(1-(1/(2^x)), add = TRUE, 
      col = "deepskyblue3",
      lwd = 2)
points (H[10]~gens[10])
points (P_Q[10]~gens[10])

legend(x = "right",             # Position
       legend = c("H", "P+Q"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c("FireBrick", "deepskyblue3"),           # Line colors
       lwd = 2)                 # Line width


### multiple loci

# A, B, C
curve(((1-(1/(2^x)))*(1-(1/(2^x)))*(1-(1/(2^x)))), from=0, to=5, main = "Genotype Proportions per \nGeneration of Inbreeding", 
      xlab="Generation (Fx)", 
      ylab="P total",
      col = "darkolivegreen",
      lwd = 2)
curve(1-(1/(2^x)), add = TRUE, 
      col = "deepskyblue3",
      lwd = 2)
legend(x = "bottomright",             # Position
       legend = c("P_total", "P_A"),  # Legend texts
       lty = c(1, 1),           # Line types
       col = c("darkolivegreen", "deepskyblue3"),           # Line colors
       lwd = 2) 


# variable number of causal loci
curve(((1-(1/(2^x)))*(1-(1/(2^x)))*(1-(1/(2^x)))), from=0, to=20, main = "Genotype Proportions per \nGeneration of Inbreeding", 
      xlab="Generation (Fx)", 
      ylab="P total",
      col = "darkolivegreen",
      lwd = 2)
for(i in 1:100){
  curve(((1-(1/(2^x)))^i), add = TRUE, 
      col = i,
      lwd = 2)
}




## Histogram of phenotype data

pericarp <- read.csv("pericarp_trimmed.csv")

# pericarp <- drop_na(pericarp)
pericarp$mean14 <- as.double(pericarp$mean14)

hist(pericarp$mean14, breaks = 20, main = "Periparp Hardness", xlab = "Force to Penetrate (N*100)", col = "darkseagreen")



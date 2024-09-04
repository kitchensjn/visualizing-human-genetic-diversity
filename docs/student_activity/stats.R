library(ggplot2)

guesses <- read.csv("guesses.csv")


ggplot(data=guesses) +
  geom_vline(xintercept=46, color="red", linewidth=2) +
  geom_histogram(aes(overlap), bins=20) +
  theme_minimal() +
  theme(text = element_text(size=30))

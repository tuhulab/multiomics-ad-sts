library(readxl)
library(dplyr)
library(ggplot2)

# examine sample sorting -----------------------------------------------
date <- 
  readr::read_csv("../multiomics-ad-transcriptomics/data/RNAseq_sample_annotation(extensive).csv") %>% 
  select(subject, visit, date_visit, skin_type) %>% 
  mutate(date_visit = date_visit %>% as.Date(format = "%d.%m.%Y")) %>%
  distinct() %>% 
  mutate(date_visit = date_visit - (date_visit %>% min))

sts_sorting <- readxl::read_excel("data/stp_sorting.xlsx") %>% 
  filter(is.na(sample_id) == FALSE) %>% 
  mutate(visit = sample_id %>% stringr::str_extract("^\\d{2}"),
         group = sample_id %>% stringr::str_extract("AD|CO"),
         subject = sample_id %>% stringr::str_extract("(AD|CO)\\_\\d{2}"),
         skin_type = sample_id %>% stringr::str_extract("LS|NL")) %>% 
  mutate(skin_type = ifelse(group == "CO", "NN", skin_type)) %>% 
  left_join(date, by = c("visit", "subject", "skin_type"))

sts_sorting %>% arrange(subject, visit, skin_type) %>% View 
# missing_sample <-
#   readr::read_csv("../multiomics-ad-transcriptomics/data/RNAseq_sample_annotation(extensive).csv") %>% 
#   select(subject, visit, date_visit, skin_type) %>% 
#   distinct() %>% 
#   mutate(date_visit = date_visit %>% as.Date(format = "%d.%m.%Y")) %>% 
#   left_join(sts_sorting, by = c("visit", "subject", "skin_type")) %>% 
#   filter(is.na(box) == TRUE) %>% select(subject, visit, date_visit = date_visit.x) %>% 
#   distinct() %>% arrange(subject, visit)
# 
# missing_sample %>% xlsx::write.xlsx("data/tape_striping_missing_samples.xlsx")

# potential batch effect
sts_sorting %>% 
  mutate(loc = factor(paste(box, column, sep = "_"), 
                      levels = paste(box, column, sep = "_"))) %>% 
  ggplot(aes(loc, date_visit)) + geom_point(aes(colour = factor(box),
                                                shape = skin_type)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = .5)) +
  scale_x_discrete(breaks = paste(1:23, "E", sep="_"),
                   labels = paste0("Box", 1:23)) +
  xlab("Location") + ylab("Date of Visit") + 
  ggtitle("Potential Batch Effect")



df2a <- df2 %>% dplyr::mutate(month = as.numeric(month),
                              month1 = dplyr::case_when(month %in% c(4, 5) ~ "Spring",
                                                        month %in% c(6, 7, 8) ~ "Summer",
                                                        TRUE ~ "Autumn")) 


p2 <- df2a %>% 
  dplyr::filter(sex != "Unknown",
                species_name == "Pipistrellus nathusii") %>% 
  ggplot(., aes(x = factor(month1, level = c("Spring", "Summer", "Autumn")),
                y = bat_capture_rate_per_trap_per_hour, colour = sex, group = sex))+
  stat_summary(geom = "point")+
  stat_smooth(method = "loess", alpha = 0.4)+
  scale_colour_manual(values = c("#332288", "#DDCC77"),
                      breaks = c("Female", "Male"))+
  theme_ipsum(base_size = 11,
              axis_title_size = 11)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.5, "lines"))+
  labs(x = "Season",
       y = "Bat capture rate (bats per trap, per hour) Mean ± 95% CI",
       colour = "Sex")

p2_build <- ggplot_build(p2)

df2 %>% 
  dplyr::filter(sex != "Unknown",
                species_name == "Pipistrellus nathusii") %>%  
  ggplot(., aes(x = month, y = bat_capture_rate_per_trap_per_hour, colour = sex))+
  stat_smooth(method = "loess", alpha = 0.4)+
  scale_x_continuous(limits = c(1, 12),
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
                     label = c("January", "February", "March", "April", "May", "June",
                               "July", "August", "September", "October", "November", "December"))+
  scale_colour_manual(values = c("#332288", "#DDCC77"),
                      breaks = c("Female", "Male"))+
  theme_ipsum(base_size = 11,
              axis_title_size = 11)+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.5, "lines"))+
  labs(x = "Month",
       y = "Bat capture rate (bats per trap, per hour) Mean ± 95% CI",
       colour = "Sex")+
  geom_hline(aes(yintercept = 0.035247778), colour = "green")+
  geom_hline(aes(yintercept = 0.018556269), colour = "red")+
  geom_hline(aes(yintercept = 0.077798612), colour = "orange")

# F:M ratio

0.035247778/ 0.146097857 # spring
0.018556269/ 0.130310516 # summer
0.077798612/ 0.211434943 # autumn

# (mean autumn female-mean summer female)/(mean autumn male-mean summer male)
(0.077798612 - 0.018556269) / (0.211434943 - 0.130310516)













gam1 <- gam::gam(bat_capture_rate_per_trap_per_hour ~ s(month, df = 6) + sex ,data = df2a)
summary(gam1)

par(mfrow = c(1, 3)) #to partition the Plotting Window
plot(gam1, se = TRUE) 
#se stands for standard error Bands
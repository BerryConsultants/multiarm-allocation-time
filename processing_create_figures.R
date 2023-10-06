#########################################################################################################
# processing_create_figures.R: Process simulation results to create manuscript figures
# Copyright (C) 2023 Berry Consultants
# 
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU 
# General Public License as published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with this program. 
# If not, see <https://www.gnu.org/licenses/>.
#########################################################################################################

library(dplyr)
options(dplyr.summarise.inform = FALSE)
library(ggplot2)
library(colorspace)

sim_folder = "SimOutput"
fig_folder = "Figures"

# Set up colors for figures
cols = qualitative_hcl(n = 5, palette = "Dark 3")[c(1,3,4,2,5)]
cols2 = qualitative_hcl(n = 5, palette = "Dark 3")[c(1,3,4,2,2,5)]
barcols = c("black", "seagreen", "mediumpurple1", 'mediumpurple4', 'steelblue', 'dodgerblue4', 'firebrick')

# Read in output file
out = readRDS(paste0(sim_folder, "/OutDF.RDS"))

out = out %>%
  filter(spendingParam == -0.5) %>% 
  mutate(DesignLabel = case_when(Design == "armdroppbo" & fixedRatio == "yes" ~ paste0("AD PBO B"), 
                                 Design == "armdroppbo" & fixedRatio == "no" ~ paste0("AD PBO A"), 
                                 Design == "armdropmax" & fixedRatio == "yes" ~ paste0("AD MAX B"),
                                 Design == "armdropmax" & fixedRatio == "no" ~ paste0("AD MAX A"), 
                                 Design == "rar" ~ "RAR", 
                                 Design == "mams" ~ "MAMS", 
                                 Design == "fixed" ~ "Fixed"))

out$VSR = factor(out$VSR, levels = c("null", 'threemixed', 'halfnugget', 'twolow', 'mixed', 'twohigh', 'nugget'))
out$VSRLabel = factor(out$VSR, labels = c("Null", "Three Mixed", "Least Favorable", 
                                          "Two Low", "Mixed", "Two High", "Nugget"))
out$TimeLabel = factor(out$TimeTrend, levels = c("flat", "linearUp", "linearDown", "seasonal", "changepoint"),
                       labels = c("Flat", "Linear Up", "Linear Down", "Seasonal", "Changepoint"))
out$fitLabel = factor(out$fitTime, levels = c("no", 'yes'), labels = c("No Time Adjustment", "Time Adjustment"))


# Create long version of output file for plotting

out.m = out %>% tidyr::pivot_longer(-c(VSR:fixedRatio, alpha, OriginalLocation, ends_with("Label")), 
                                    names_to = "Metric", values_to = "value")
out.m$Design = paste0(out.m$Design, ifelse(out.m$activetowin, "_ATW", ""))
out.m$Design = factor(out.m$Design, levels = c("rar", "detrar", "armdropmax", 
                                               "armdropmax_ATW", "armdroppbo", 
                                               "armdroppbo_ATW", "mams_ATW", "fixed"))

out.m = out.m %>% filter(
  Metric %in% c("power", "effectMSE", "regret", 'avgSSSelectedArm', "idp", 'bestArmRate'),
  Design %in% c("rar", "fixed", "armdropmax_ATW", "armdropmax", "armdroppbo_ATW", "mams_ATW")) %>%
  mutate(MetricLabel = factor(Metric, 
                              levels = c("power", "avgSSSelectedArm", "idp", "bestArmRate", 
                                         "expectedResponders", "effectMSE", "regret"),
                              labels = c("Power", "Avg Sample Size Selected Arm","IDP",  
                                         "Best Arm Selected", "Expected Responders", "Effect MSE", 
                                         "Responders Below Optimal (Regret)")))

# Create blank_data dataframe to control the y-axis limits w/ facet_wrap
lim.power.nonnull = c(min(out$power[out$VSR!="null"]), max(out$power[out$VSR!="null"]))

tmp = out.m %>% select(-value) %>% unique()
min_data = tmp %>%
  mutate(value = case_when(Metric == "power" ~ min(out$power[out$VSR != "null"], na.rm = TRUE),
                           Metric == "idp" ~ min(out$idp[out$VSR != "null"], na.rm = TRUE),
                           Metric == "regret" ~ min(out$regret[out$VSR != "null"], na.rm = TRUE),
                           Metric == "effectMSE" ~ min(out$effectMSE[out$effectMSE != "null"], na.rm = TRUE), 
                           Metric == "expectedResponders" ~ min(out$expectedResponders[out$VSR != "null"], na.rm = TRUE), 
                           Metric == "bestArmRate" ~ min(out$bestArmRate[out$VSR != "null"], na.rm = TRUE), 
                           Metric == "avgSSSelectedArm" ~ min(out$avgSSSelectedArm[out$VSR != "null"], na.rm = TRUE)))
max_data = tmp %>%
  mutate(value = case_when(Metric == "power" ~ 1.1*max(out$power[out$VSR != "null"], na.rm = TRUE),
                           Metric == "idp" ~ 1.1*max(out$idp[out$VSR != "null"], na.rm = TRUE),
                           Metric == "regret" ~ 1.1*max(out$regret[out$VSR != "null"], na.rm = TRUE),
                           Metric == "effectMSE" ~ 1.1*max(out$effectMSE[out$VSR != "null"], na.rm = TRUE), 
                           Metric == "expectedResponders" ~ 1.1*max(out$expectedResponders[out$VSR != "null"], na.rm = TRUE), 
                           Metric == "bestArmRate" ~ 1.1*max(out$bestArmRate[out$VSR != "null"], na.rm = TRUE), 
                           Metric == "avgSSSelectedArm" ~ 1.1*max(out$avgSSSelectedArm[out$VSR != "null"], na.rm = TRUE)))
blank_data = bind_rows(min_data, max_data)

out.m <- out.m %>%
  mutate(value_rounded = case_when(Metric %in% c("power", "effectMSE", "bestArmRate") ~ 
                                     sprintf("%0.0f", 100*value), 
                                   TRUE ~ sprintf("%0.0f", value)), 
         value_rounded2 = case_when(Metric == "power" ~ sprintf("%0.1f", value*100)), 
         DesignLabel = factor(DesignLabel, levels = c("Fixed", "MAMS", "AD PBO A", "AD PBO B", 
                                                      "AD MAX A", "AD MAX B", "RAR")))

###############################################
## Create plots for manuscript 
vsr_unique = unique(out$VSR[out$VSR!="null"])
vsr_labels = unique(out$VSRLabel[out$VSR!="null"])
fit_time = c('yes', 'no')
time_labels = c("Time Adjustment", "No Time Adjustment")

# OC figures when no time adjustment in model
pdf(file = paste0(fig_folder, '/OC_Figures_fitTimeNo.pdf'), 10, 10)
for(v in 1:length(vsr_unique)){
  for(t in 2){
    title = paste0(vsr_labels[v], " Efficacy Scenario; ", time_labels[t])
    
    print(ggplot(out.m %>% filter(VSR == vsr_unique[v], fitTime == fit_time[t]), 
                 aes(x = TimeLabel, fill = DesignLabel, color = DesignLabel)) +
            geom_text(aes(y = value*1.05, label = value_rounded), 
                      size = 2, position = position_dodge(0.85)) +
            geom_bar(aes(y = value), stat = "identity", 
                     position = position_dodge(0.85), width = 0.75) + 
            scale_color_manual(values = barcols) +
            scale_fill_manual(values = barcols) +
            guides(color = 'none') + 
            facet_wrap(~MetricLabel, scales = 'free', ncol = 2) + 
            theme_bw() + 
            labs(x = "Time Scenario", 
                 title = title, 
                 y = "", 
                 fill = "") + 
            geom_blank(aes(y = value), data = blank_data) + 
            theme(axis.title.x = element_text(size = 11), 
                  strip.text = element_text(size = 12), 
                  legend.position = "bottom"))
    
  }
}
dev.off()

# OC figures when time adjustment in model
pdf(file = paste0(fig_folder, '/OC_Figures_fitTimeYes.pdf'), 10, 10)
for(v in 1:length(vsr_unique)){
  for(t in 1){
    title = paste0(vsr_labels[v], " Efficacy Scenario; ", time_labels[t])
    
    print(ggplot(out.m %>% filter(VSR == vsr_unique[v], fitTime == fit_time[t]), 
                 aes(x = TimeLabel, fill = DesignLabel, color = DesignLabel)) +
            geom_text(aes(y = value*1.05, label = value_rounded), 
                      size = 2, position = position_dodge(0.85)) +
            geom_bar(aes(y = value), stat = "identity", 
                     position = position_dodge(0.85), width = 0.75) + 
            scale_color_manual(values = barcols) +
            scale_fill_manual(values = barcols) +
            guides(color = 'none') + 
            facet_wrap(~MetricLabel, scales = 'free', ncol = 2) + 
            theme_bw() + 
            labs(x = "Time Scenario", 
                 title = title, 
                 y = "", 
                 fill = "") + 
            geom_blank(aes(y = value), data = blank_data) + 
            theme(axis.title.x = element_text(size = 11), 
                  strip.text = element_text(size = 12), 
                  legend.position = "bottom"))
    
  }
}
dev.off()


pdf(file = paste0(fig_folder, '/OC_Figures_fitTimeNo_typeIerror.pdf'), 6, 4)

title = paste0("Null Efficacy Scenario; ", time_labels[2])
print(ggplot(out.m %>% filter(VSR == "null", MetricLabel == "Power", fitTime == fit_time[2]), 
             aes(x = TimeLabel, fill = DesignLabel, color = DesignLabel)) +
        geom_text(aes(y = value+0.003, label = value_rounded2), 
                  size = 2.4, position = position_dodge(0.85)) +
        geom_bar(aes(y = value), stat = "identity", 
                 position = position_dodge(0.85), width = 0.75) + 
        scale_color_manual(values = barcols) +
        scale_fill_manual(values = barcols) +
        guides(color = 'none') + 
        theme_bw() + 
        labs(x = "Time Scenario", 
             title = title, 
             fill = "") + 
        scale_y_continuous("Simulated type I error rate", breaks = seq(0, 0.1, 0.025), limits = c(0, 0.1)) +
        theme(axis.title.x = element_text(size = 11), 
              strip.text = element_text(size = 12), 
              legend.position = "bottom"))
dev.off()

pdf(file = paste0(fig_folder, '/OC_Figures_fitTimeYes_typeIerror.pdf'), 6, 4)
title = paste0("Null Efficacy Scenario; ", time_labels[1])
print(ggplot(out.m %>% filter(VSR == "null", MetricLabel == "Power", fitTime == fit_time[1]), 
             aes(x = TimeLabel, fill = DesignLabel, color = DesignLabel)) +
        geom_text(aes(y = value+0.003, label = value_rounded2), 
                  size = 2.4, position = position_dodge(0.85)) +
        geom_bar(aes(y = value), stat = "identity", 
                 position = position_dodge(0.85), width = 0.75) + 
        scale_color_manual(values = barcols) +
        scale_fill_manual(values = barcols) +
        guides(color = 'none') + 
        # facet_wrap(~MetricLabel, scales = 'free', ncol = 2) + 
        theme_bw() + 
        labs(x = "Time Scenario", 
             title = title, 
             fill = "") + 
        # geom_blank(data = blank_data) + 
        scale_y_continuous("Simulated type I error rate", breaks = seq(0, 0.1, 0.025), limits = c(0, 0.1)) +
        theme(axis.title.x = element_text(size = 11), 
              strip.text = element_text(size = 12), 
              legend.position = "bottom"))
dev.off()


pdf(file = paste0(fig_folder, 'OC_Figures_fitTimeNo_timeScenarioFlat.pdf'), 10, 10)

print(ggplot(out.m %>% filter(VSR != "null", fitTime == 'no', TimeTrend == "flat"), 
             aes(x = VSRLabel, fill = DesignLabel, color = DesignLabel)) +
        geom_text(aes(y = value*1.05, label = value_rounded), 
                  size = 2, position = position_dodge(0.85)) +
        geom_bar(aes(y = value), stat = "identity", 
                 position = position_dodge(0.85), width = 0.75) + 
        scale_color_manual(values = barcols) +
        scale_fill_manual(values = barcols) +
        guides(color = 'none') + 
        facet_wrap(~MetricLabel, scales = 'free', ncol = 2) + 
        theme_bw() + 
        labs(x = "Efficacy Scenario", 
             title = "Flat Time Scenario",  
             y = "", 
             fill = "") + 
        geom_blank(data = blank_data %>% filter(VSR != "null")) +
        theme(axis.text.x = element_text(angle = 15, hjust=1, size = 11), 
              axis.title.x = element_text(size = 11), 
              strip.text = element_text(size = 12), 
              legend.position = "bottom"))
dev.off()

pdf(file = paste0(fig_folder, '/OC_Figures_null_noTimeAdj.pdf'), 10, 5)
print(ggplot(out.m %>% filter(VSR == "null", fitTime == 'no', 
                              Metric %in% c("avgSSSelectedArm", "effectMSE")), 
             aes(x = TimeLabel, fill = DesignLabel, color = DesignLabel)) +
        geom_text(aes(y = value*1.05, label = value_rounded), 
                  size = 2, position = position_dodge(0.85)) +
        geom_bar(aes(y = value), stat = "identity", 
                 position = position_dodge(0.85), width = 0.75) + 
        scale_color_manual(values = barcols) +
        scale_fill_manual(values = barcols) +
        guides(color = 'none') + 
        facet_wrap(~MetricLabel, scales = 'free', ncol = 2) + 
        theme_bw() + 
        labs(x = "Time Scenario", 
             title = "Null Efficacy Scenario; No Time Adjustment",  
             y = "", 
             fill = "") + 
        theme(axis.text.x = element_text(angle = 15, hjust=1, size = 11), 
              axis.title.x = element_text(size = 11), 
              strip.text = element_text(size = 12), 
              legend.position = "bottom"))
dev.off()

tmp <- out.m %>% filter(VSR != "null", TimeTrend == "flat")
df.time <- tmp %>%
  group_by(DesignLabel, VSRLabel, MetricLabel, TimeTrend) %>%
  summarize(metric_change = (value[fitTime == "yes"] - value[fitTime == "no"])) %>%
  mutate(TimeLabel = factor(case_when(TimeTrend == "flat" ~ "Flat",
                                      TimeTrend == "linearUp" ~ "Linear Up",
                                      TimeTrend == "linearDown" ~ "Linear Down",
                                      TimeTrend == "seasonal" ~ "Seasonal",
                                      TimeTrend == "changepoint" ~ "Changepoint"),
                            levels = c("Flat",
                                       "Linear Up",
                                       "Linear Down",
                                       "Seasonal",
                                       "Changepoint"))) %>%
  mutate(DesignLabel = factor(DesignLabel, levels = c("Fixed", "MAMS", "AD PBO A", "AD PBO B", 
                                                      "AD MAX A", "AD MAX B", "RAR")))
mindf <- df.time %>%
  group_by(DesignLabel, VSRLabel, MetricLabel) %>%
  summarize(metric_change = -1*max(abs(metric_change)))
maxdf <- df.time %>%
  group_by(DesignLabel, VSRLabel, MetricLabel) %>%
  summarize(metric_change = max(abs(metric_change)))
blank_dat <- bind_rows(mindf, maxdf)

p1 <- (ggplot(df.time,
              aes(x = VSRLabel,
                  y = metric_change,
                  fill = DesignLabel)) +
         geom_bar(aes(y = metric_change), stat = "identity", 
                  position = position_dodge(0.85), width = 0.75) + 
         geom_hline(yintercept = 0, lty = 2) +
         scale_color_manual(values = barcols) +
         scale_fill_manual(values = barcols) +
         guides(color = 'none') +
         facet_wrap(~MetricLabel, ncol = 2, scales = 'free_y') +
         theme_bw() +
         labs(x = "Efficacy Scenario",
              title = "Flat Time Scenario",
              y = "Change in metric (Time covariate versus unadjusted)",
              fill = "Design") +
         geom_blank(data = blank_dat) +
         theme(legend.position = 'bottom',
               axis.text.x = element_text(angle = 15, hjust=1, size = 10),
               strip.text = element_text(size = 12)) +
         theme(
           strip.text.y = element_blank()
         ))

pdf(file = paste0(fig_folder, '/OC_Figures_ChangeInMetrics.pdf'), 10, 10)
p1
dev.off()


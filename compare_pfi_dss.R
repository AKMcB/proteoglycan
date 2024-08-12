
dss <- as.data.frame(fread("cofficient_sign_dss.csv"))
pfi <- as.data.frame(fread("cofficient_sign_pfi.csv"))

dss <- dss[,-1]
pfi <- pfi[,-1]

dss$analysis <- "DSS"
pfi$analysis <- "PFI"

dss <- subset(dss, dss$coef_sign == ">0")
pfi <- subset(pfi, pfi$coef_sign == ">0")

combined<- rbind(dss,pfi)




wide_df <- pivot_wider(data = combined,
                       names_from =coef_sign,
                       values_from = exp.coef.)



attach(wide_df)
# Define size breaks and labels
size_breaks <- c(0.5,1,2,4)
size_labels <- as.character(size_breaks)

p <- ggplot(data= wide_df, aes(x=analysis, y=Gene))+ 
  geom_point(aes(size=`>0`, color= analysis)) +
  theme_bw()+
  scale_color_manual(values = c("DSS" = "orange", "PFI" = "darkblue"))+
  scale_size_continuous(
    range = c(0.5, 4),  # Adjust the range as needed
    breaks = size_breaks,
    labels = size_labels)+ 
  labs(x = "Survival Analysis", y = "Genes", 
       size = "Hazard Ratio",colour = "Survival Analysis")+
  ggtitle("Increased Hazard: DSS & PFI")+
  facet_wrap(~Cancer, ncol = 10)
p


pdf("dotplot_survival_pfi_dss_increased_hazard.pdf", width = 8, height = 15, onefile = F)
print(p)
dev.off()

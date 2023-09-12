#######################################################
# Define colour palette (CMS_pal)
## Carolin Sauer
#######################################################


# Colour Codes
CMSpink <- "#B74271FF"
CMSyellow <- "#E6B800FF"
CMSgrey <- "#838377FF"
CMSteal <- "#008080FF"
CMSpurple <- "#4E2574FF"
CMSblue <- "#1C3C6CFF"
CMSlightpink <- "#E287B3FF"
CMSgreen <- "#236050FF"
CMSlightblue <- "#71A7B3FF"
CMSred <- "#840013FF"
CMSorange <- "#cc5200"
CMSlightgreen <- "#339933"
CMSlightred <- "#b30000"
CMSdarkgrey <- "#333333"
CMSroyalblue <- "#1c41b0"


# original colour palette
CMS_pal <- c(CMSblue,CMSgreen,CMSgrey,CMSlightblue,CMSlightpink,CMSpink,CMSpurple,CMSred,CMSteal,CMSyellow, CMSorange, CMSlightgreen, CMSlightred, CMSdarkgrey, CMSroyalblue)

show_CMS_pal <- function(palette){
  barplot(c(1:length(palette)), col = palette)
}





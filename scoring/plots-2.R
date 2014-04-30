require(data.table)
require(ggplot2)


read_data <- function(fname) {
    d <- fread(fname)
    setnames(d, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"),
                c("aa", "ca", "cb", "vcb", "cos_theta", "theta", "d_theta", "cos_omega", "omega", "d_omega"))
    setkey(d, "aa")
    d[omega == max(d$omega)]$omega <- 0
    d$omega <- 2*d$omega
    return(d)
}

aa_list <- c("ALA", "LEU", "PRO", "GLY", "ASP", "ASN", "TYR", "HIS", "GLU", "CYS", "PHE", "VAL", "ILE", "ARG", "THR", "LYS", "SER", "GLN", "MET", "TRP")

nat <- read_data("distance_stats_2.tab")
back <- read_data("distance_stats_2_shuffled.tab")
pairs <- levels(as.factor(nat$aa))

from =  2.0
to   = 30.0
from_show = 2.5
to_show = 25.0

plot_density <- function(aa, data, xparam, label) {

    filename <- paste("plots/potential/", xparam, "/", label, "/", aa, ".png", sep="")
    # if (is.na(data[aa][1, tolower(xparam), with=FALSE])) {
    #     print(paste("No data for", aa))
    #     file.remove(filename)
    #     return()
    # }

    print(aa)
    d <- density(data[aa][[tolower(xparam)]], from=from, to=to)
    dens <- as.data.frame(cbind(
        x= d$x,
        density = d$y))

    p <- ggplot(aes(x=x, y=density), data=dens) +
      geom_line() +
      # xlim(from_show, to_show) +
      labs(title=paste("Probability density ", aa), x=xparam) +
      theme_classic()
    ggsave(filename, width=5, height=5)
}


plot_density_ratio <- function(aa, nat, back, xparam, label) {
    filename <- paste("plots/potential/", xparam ,"/", label, "/", aa, ".png", sep="")
    # if (is.na(nat[aa][1]$ca) || is.na(back[aa][1]$ca)) {
    #     print(paste("No data for", aa))
    #     file.remove(filename)
    #     return()
    # }

    print(aa)

    d.nat <- density(nat[aa][[tolower(xparam)]], from=from, to=to)
    d.back <- density(back[aa][[tolower(xparam)]], from=from, to=to)

    d.logratio <- as.data.frame(cbind(
        x = d.nat$x,
        logodds = log(d.nat$y / d.back$y)))

    p <- ggplot(aes(x=x, y=logodds), data = d.logratio) +
      geom_line() +
      geom_hline(yintercept=0, color="red") +
      # xlim(from_show, to_show) +
      # ylim(from_show, to_show) +
      labs(title=paste("Log odds", aa), x=xparam) +
      theme_classic()
    ggsave(filename, width=5, height=5)
}


############### 
plot_all <- function(xparam) {
    dev.null <- lapply(pairs, plot_density, nat, xparam, "nat")
    dev.null <- lapply(pairs, plot_density, back, xparam, "back")

    dev.null <- lapply(pairs, plot_density_ratio, nat, back, xparam, "logratio")
}

dev.null <- lapply(c("Cb", "vCb", "Ca"), plot_all)


# > plot(d["TYR_TYR"]$theta, d["TYR_TYR"]$dCb)
# > d.close <- d[dvCb >0][dvCb < 10]
# > hist(d.close["ASP_ASP"]$dvCb)
# > quartz()
# > hist(d.close["ASP_LYS"]$dvCb)
# > hist(d.close["ASP_ARG"]$dvCb)
# > hist(d.close["ARG_ARG"]$dvCb)
# > hist(d.close["ARG_LYS"]$dvCb)
# > hist(d.close["HIS_HIS"]$dvCb)
# > hist(d.close["TYR_TYR"]$dvCb)
# > hist(d.close["CYS_CYS"]$dvCb)
# > hist(d.close["GLU_GLU"]$dvCb)
# > hist(d.close["GLU_ARG"]$dvCb)
# > hist(d.close["GLU_GLU"]$dvCb)
# > hist(d.close["GLU_ARG"]$dvCb)
# > hist(d.close["GLU_GLU"]$dvCb)
# > hist(d.close["GLU_ARG"]$dvCb)
# > hist(d.close["GLU_GLU"]$dvCb)
# > hist(d.close["SER_SER"]$dvCb)
# > hist(d.close["THR_THR"]$dvCb)
# > hist(d.close["GLY_GLY"]$dvCb)
# Error in hist.default(d.close["GLY_GLY"]$dvCb) :
#   invalid number of 'breaks'
# > hist(d.close["ILE_ILE"]$dvCb)
# > hist(d.close["VAL_VAL"]$dvCb)
# > hist(d.close["PHE_PHE"]$dvCb)
# > hist(d.close[7<dvCb<8]$dvCb)
# Error: unexpected '<' in "hist(d.close[7<dvCb<"
# > hist(d.close[dvCb<8 && dvCb >7]$dvCb)
# Error in hist.default(d.close[dvCb < 8 && dvCb > 7]$dvCb) :
#   invalid number of 'breaks'
# > hist(d.close[dvCb<8 & dvCb >7]$dvCb)
# > length(d.close[dvCb<8 & dvCb >7]$dvCb)
# [1] 357024
# > length(d.close[dvCb<8 & dvCb >7]$dvCb) / length(d.close$dvCb)
# [1] 0.1700473
# > length(d.close[ASP_ASP][dvCb<8 & dvCb >7]$dvCb) / length(d.close[ASP_ASP]$dvCb)
# Error in eval(expr, envir, enclos) : object 'ASP_ASP' not found
# > length(d.close["ASP_ASP"][dvCb<8 & dvCb >7]$dvCb) / length(d.close["ASP_ASP"]$dvCb)
# [1] 0.1681183
# > length(d.close[dvCb>4 & dvCb <5]$dvCb) / length(d.close$dvCb)
# [1] 0.05926584
# > length(d.close["ASP_ASP"][dvCb>4 & dvCb <5]$dvCb) / length(d.close["ASP_ASP"]$dvCb)
# [1] 0.05544034
# > length(d.close["ASP_ARG"][dvCb>4 & dvCb <5]$dvCb) / length(d.close["ASP_ARG"]$dvCb)


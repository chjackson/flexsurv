## Code to generate the hex logo

library(ggplot2)
library(flexsurv)
library(hexSticker)
#remotes::install_github("AllanCameron/geomtextpath")
sysfonts::font_add(family = "Aller_Rg", regular = "Aller_rg.ttf")
showtext::showtext_auto()

dat <- data.frame(
    x = c(0, 1, 5, 12, 15),
    y = c(1, 0.7, 0.6, 0.2, 0.2)
)
sdat <- data.frame(
    x = c(0, 2,   4,   6,     10,   12,  12.5,  15),
    y = c(1, 0.7, 0.65, 0.62, 0.55, 0.50, 0.25, 0.15)
)
sfit <- smooth.spline(sdat$x, sdat$y, spar=0.2)
sfit <- as.data.frame(predict(sfit, x = seq(0,15,by=0.1)))

textdf <- data.frame(
    x = c(0,   2,   4,   6,     10,   12,  12.5,  15),
    y = c(0.82, 0.68, 0.65, 0.62, 0.55, 0.50, 0.4, 0.35)
)
textdf <- smooth.spline(textdf$x, textdf$y, spar=0.2)
textdf <- as.data.frame(predict(textdf))
textdf$text <- "flexsurv"

fsp <- ggplot(dat, aes(x=x,y=y)) + 
    xlab("") + ylab("") +
    theme(plot.background = element_rect(fill="transparent", colour="transparent"), 
          panel.background = element_rect(fill="transparent"), 
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.margin = unit(c(0, 0, -0.05, -0.05), "null"),
          panel.margin = unit(c(0, 0, -0.05, -0.05), "null"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())+
    geom_step(size=2, col="gray60") + 
    geom_line(data=sfit, col="darkblue", size=2.3, alpha=0.7) + 
 #   geom_segment(aes(y = 0, x=0, yend=0, xend=15), col="gray70") +
 #   geom_segment(aes(x=0, xend=0, y=0, yend=1.1), col="gray70") + 
    geomtextpath::geom_textpath(data=textdf, aes(label=text), 
                                size=7, vjust=-0.5, hjust=0, 
                                text_only=TRUE,
                                col="gray20",
                                family = "Aller_Rg", 
                                fontface="bold")  +
    coord_cartesian(ylim=c(0,1), xlim=c(0,15), clip="off") 
fsp

## from https://github.com/GuangchuangYu/hexSticker

theme_sticker <- function(size=1.2, ...) {
    center <- 1
    radius <- 1
    h <- radius
    w <- sqrt(3)/2 * radius
    m <- 1.02
    list(
        theme_transparent() +
            theme(plot.margin = margin(b = -.2, l= -.2, unit = "lines"),
                  strip.text = element_blank(),
                  line = element_blank(),
                  text = element_blank(),
                  title = element_blank(), ...),
        coord_fixed(),
        scale_y_continuous(expand = c(0, 0), limits = c(center-h*m , center+h*m )),
        scale_x_continuous(expand = c(0, 0), limits = c(center-w*m , center+w*m ))
    )
}

p <- ggplot() + 
    geom_hexagon(size = 2, fill = "azure", color = "gray20") + 
    geom_subview(subview=fsp, x=0.9, y=0.8, width=1.5, height=1.5) +
    theme_sticker()
p

showtext::showtext_opts(dpi = 300)
ggsave(p, width=43.9, height=50.8, units="mm", device=png,
       bg="transparent", filename="man/figures/flexsurv_hex.png")


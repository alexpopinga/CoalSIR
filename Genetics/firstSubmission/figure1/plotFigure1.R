require(rjson)

df <- fromJSON(file='SIR_output.json')

Scol <- rgb(0,0,0)
Icol <- rgb(0.9,0.4,0.4)
Rcol <- rgb(0,0.8,0.8)

pdf('figure1.pdf', width=7, height=6)
plot(df$t, df$S, 's', lwd=2, col=Scol, ylim=c(0,1000),
     xlab='Time', ylab='Population size')
lines(df$t, df$I, 's', lwd=2, col=Icol)
lines(df$t, df$R, 's', lwd=2, col=Rcol)

legend('left', inset=.05, c('S(t)','I(t)','R(t)'), lwd=2,
       col=c(Scol, Icol, Rcol))

dev.off()

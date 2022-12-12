### NO AGN ### 

plot.ProSpectSEDlike=function(x, xlim=c(1e3,1e7), ylim='auto',
                              xlab='Wavelength (Ang)', ylab='auto', grid=TRUE,
                              type='flux', ...){
  
  .cgs_to_jansky = 1e23

  StarFlux=Lum2Flux(wave=x$SEDout$Stars$wave_lum, lum=x$SEDout$Stars$lum_unatten, z=x$Data$arglist$z,  ref= 'Planck15', LumDist_Mpc=x$Data$arglist$LumDist_Mpc)
  StarFlux$flux=convert_wave2freq(flux_wave=StarFlux$flux*.cgs_to_jansky, wave=StarFlux$wave)
  
  AttenStarFlux=Lum2Flux(wave=x$SEDout$Stars$wave_lum, lum=x$SEDout$Stars$lum_atten, z=x$Data$arglist$z,  ref= 'Planck15', LumDist_Mpc=x$Data$arglist$LumDist_Mpc)
  AttenStarFlux$flux=convert_wave2freq(flux_wave=AttenStarFlux$flux*.cgs_to_jansky, wave=StarFlux$wave)
  
  DustFlux=Lum2Flux(wave=x$SEDout$DustEmit$wave, lum=x$SEDout$DustEmit$lum, z=x$Data$arglist$z,  ref= 'Planck15', LumDist_Mpc=x$Data$arglist$LumDist_Mpc)
  DustFlux$flux=convert_wave2freq(flux_wave=DustFlux$flux*.cgs_to_jansky, wave=x$SEDout$DustEmit$wave)
  
  if(!is.null(x$SEDout$AGN)){
    AGNFlux=Lum2Flux(wave=x$SEDout$AGN$wave, lum=x$SEDout$AGN$lum, z=x$Data$arglist$z,  ref= 'Planck15', LumDist_Mpc=x$Data$arglist$LumDist_Mpc)
    AGNFlux$flux=convert_wave2freq(flux_wave=AGNFlux$flux*.cgs_to_jansky, wave=AGNFlux$wave)
  
    FullFlux = approxfun(x$SEDout$FinalFlux)
    AGNFluxApprox = approxfun(AGNFlux)
  
    minwav = 5e4 * (1+x$Data$arglist$z)
    maxwav= 20e4 * (1+x$Data$arglist$z)
    AGNfrac = integral(AGNFluxApprox, xmin = minwav, xmax = maxwav)/integral(FullFlux, xmin = minwav, xmax = maxwav)
}
  
  if(type=='flux'){
    par(oma=c(3.1,3.1,1.1,1.1))
    par(mar=c(0,0,0,0))
    par(fig=c(0,0.55,0.25,1))
    plot(x$SEDout, xlim=xlim, ylim=ylim, xlab='', ylab=ylab, grid=grid, type='flux', labels = c(F,T), ...)
    points(x$Data$flux[,2:3], pch=16, col='red')
    if(requireNamespace("magicaxis", quietly=TRUE)){
      magicaxis::magerr(x$Data$flux[[2]], x$Data$flux[[3]], ylo=x$Data$flux[[4]], col='red')  
      if(!is.null(x$SEDout$AGN)){lines(AGNFlux, col = 'purple')}
      lines(StarFlux, col = 'blue')
      lines(AttenStarFlux, col = 'red')
      lines(DustFlux, col = 'brown')
    }
    
    if(!is.null(x$SEDout$AGN)){legend('bottomright', legend=c(paste('LP =',round(x$LP,3)),
                                   paste('log SMformed =', round(log10(x$SEDout$Stars$masstot), 3)),
                                   paste('SFR =', round(x$SEDout$Stars$SFRburst,3)),
                                   paste('Zfinal =', round(x$parm['Zfinal'],3)),
                                   paste('log AGNlum =', round(x$parm['AGNlum'],3)),
                                   paste('AGNfrac =', round(AGNfrac,3))))
    } else {legend('bottomright', legend=c(paste('LP =',round(x$LP,3)),
                                                                paste('log SMformed =', round(log10(x$SEDout$Stars$masstot), 3)),
                                                                paste('SFR =', round(x$SEDout$Stars$SFRburst,3)),
                                                                paste('Zfinal =', round(x$parm['Zfinal'],3))))
                                   }

    par(mar=c(0,0,0,0))
    par(fig=c(0,0.55,0,0.25), new = TRUE)
    if(requireNamespace("magicaxis", quietly=TRUE)){
      magicaxis::magplot(x$Data$flux[[2]], (x$Data$flux[[3]]-x$SEDout$Photom)/x$Data$flux[[4]], pch=16, col='red', grid=grid,
                         log='x', xlim=xlim, ylim=c(-4,4), xlab=xlab, ylab='(Data-Model)/Error')
      
    }else{
      plot(x$Data$flux[,2], x$Data$flux[,3]-x$SEDout$Photom, pch=16, col='red',
           log='x', xlim=xlim, ylim=c(-4,4), xlab=xlab, ylab='(Data-Model)/Error')
    }
    
    par(mar=c(0,0,0,0))
    par(fig=c(0.65,1,0,0.5), new = TRUE)
    if(requireNamespace("magicaxis", quietly=TRUE)){
      magicaxis::magplot(x$SEDout$Stars$agevec/1e9, x$SEDout$Stars$Zvec, type = 'l', xlab = 'Age (Gyrs)', ylab = 'Z', xlim = c(0, 13.8))  
    }
    
    par(mar=c(0,0,0,0))
    par(fig=c(0.65,1,0.5,1), new = TRUE)
    if(requireNamespace("magicaxis", quietly=TRUE)){
      magicaxis::magplot(x$SEDout$Stars$agevec/1e9, x$SEDout$Stars$SFR, type = 'l', ylab = 'SFR (Msol/yr)', xlim = c(0, 13.8), labels = c(F,T))  
    }
    
    
  }else if(type=='lum'){
    plot(x$SEDout, xlim=xlim, ylim=ylim, xlab='', ylab=ylab, grid=grid, type='lum', ...)
  }
}

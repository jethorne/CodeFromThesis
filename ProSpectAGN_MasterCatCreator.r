## Master Cat creator from highest likelihood step

# Loading all the required packages 
.libPaths(c("/group/pawsey0160/software/sles12sp3/apps/sandybridge/gcc/4.8.5/r/3.6.3/lib64/R/library",.libPaths()))
#.libPaths(c("/group/pawsey0160/jthorne/R/magnus/",.libPaths()))
library(ProSpect)
library(data.table)
library(foreach)
library(doParallel)
library(celestial)
CoreNumber = 24
registerDoParallel(cores=CoreNumber)

data(BC03lr)
data(Dale_Msol)
data(EAZY_filters)
filters = c('galex_fuv', 'galex_nuv', 'cfht_u', 'hsc_g','hsc_r','hsc_i','hsc_z', 'vista_Y', 'vista_J', 'vista_H', 'vista_Ks','irac_1', 'irac_2', 'irac_3', 'irac_4', 'mips_24', 'mips_70', 'pacs_100', 'pacs_160', 'spire_250', 'spire_350', 'spire_500')
EAZY_filter_names = c("CAPAK/galex1500.res","CAPAK/galex2500.res", "COSMOS/CFHT_filter_u.txt", "g_HSC.txt","r_HSC.txt"  , "i_HSC.txt", "z_HSC.txt" , 'VISTA/.s?_system', "IRAC/irac_tr1_2004-08-09.dat", "IRAC/irac_tr2_2004-08-09.dat", "IRAC/irac_tr3_2004-08-09.dat" ,"IRAC/irac_tr4_2004-08-09.dat", "mips/24" , "mips/70",  "herschel/pacs/100" , "herschel/pacs/160",
                      "herschel/spire/200", "herschel/spire/350", "herschel/spire/500")

EAZY_filter_locations = foreach(ii = 1:length(EAZY_filter_names), .combine = 'c')%do%{
  loc = grep(EAZY_filter_names[ii], EAZY_filters$info, ignore.case = T)
  return(loc)
}

filtout= foreach(jj = 1:length(EAZY_filter_locations))%do%{approxfun(EAZY_filters$filters[EAZY_filter_locations[jj]][[1]])}
cenwave=foreach(i =1:22, .combine = 'c')%do%{cenwavefunc(filtout[[i]])}


OutputDir = "/group/pawsey0335/jthorne/DEVILS/ProSpect/D10/Outputs/"

source('/group/pawsey0160/jthorne/DEVILS/ProSpect/D10/Code/ProSpect_Functions.R')

.lsol_to_erg= 3.828e33
.mpc_to_cm = 3.08568e+24
.lsol_to_absolute = .lsol_to_erg/(4*pi*(.mpc_to_cm*1e-5)^2)




#### Highlander Best Fit Outputs ####

foreach(jj = 1:7)%do%{ ## looping over the 7 redshift bins
  #Setting the fitting version number
  version_no = paste0('v3.42/', jj, '/')
  
  myFiles <- list.files(path= paste(OutputDir, version_no,  sep=''), pattern="*.rds")
  length(myFiles)
  
  ## looping over each galaxy output 
  Sample = foreach(i=1:length(myFiles), .combine='rbind') %dopar% {
      file_name = myFiles[i]
      filename = paste0(OutputDir,version_no, '/',file_name)
      GalStuff = readRDS(filename)
      GalID = as.character(sapply(strsplit(file_name, "\\."), "[[", 1))
      
      Redshift = GalStuff$Data$arglist$z
      
      LookbackTime = cosdistTravelTime(Redshift, ref = 'Planck15')*1e9
      agemax = 13.38e9 - LookbackTime
      
      nparm = length(GalStuff$Data$parm.names)
      post = cbind(GalStuff$Highlander_output$LD_last$Posterior1, GalStuff$Highlander_output$LD_last$Monitor)
      post = data.table(post)
      post = post[order(-LP_mon),]
      
      ## calculating a chisq cut based on the number of parameters to isolate a 1sigma range
      chisq_cut = max(post$LP_mon) - qchisq(0.68,nparm)/2
      toppost = post[post$LP_mon > chisq_cut,]
      toppost[,LP_mon:=NULL]
      
      values = as.numeric(post[1,1:nparm])
      
      BinMasses=SMstarfunc(massfunc = massfunc_snorm_trunc,
                           mSFR=10^(values[1]),
                           mpeak=(values[2]),
                           mperiod=10^(values[3]),
                           mskew=( values[4]),
                           Zfinal = 10^(values[9]),
                           magemax = agemax/1e9,
                           Z= Zfunc_massmap_lin, #yield=Data$arglist$yield,
                           z = Redshift,
                           ref = 'Planck15')
      
      ratio = as.numeric(BinMasses['TotSMstar']/toppost[1,'masstot'])
      masstotloc = which(colnames(toppost) ==  'masstot')
      rangeofvals = foreach(ii = 1:ncol(toppost), .combine = 'rbind')%do%{range(toppost[,..ii])}
      TotSMstarlimits = rangeofvals[which(colnames(toppost) == 'masstot'),] * ratio
      
      
      ## AGN fluxstuff
      x = GalStuff$ProSpectSEDlike
      .cgs_to_jansky = 1e23
      
      AGNFlux=Lum2Flux(wave=x$SEDout$AGN$wave, lum=x$SEDout$AGN$lum, z=x$Data$arglist$z,  ref= 'Planck15', LumDist_Mpc=x$Data$arglist$LumDist_Mpc)
      AGNFlux$flux=convert_wave2freq(flux_wave=AGNFlux$flux*.cgs_to_jansky, wave=AGNFlux$wave)
    
      FullFlux = approxfun(x$SEDout$FinalFlux)
      AGNFluxApprox = approxfun(AGNFlux)
      
      minwav = 5e4 * (1+x$Data$arglist$z)
      maxwav= 20e4 * (1+x$Data$arglist$z)
      AGNfrac = integral(AGNFluxApprox, xmin = minwav, xmax = maxwav)/integral(FullFlux, xmin = minwav, xmax = maxwav)
      
      ## iteration information
      type = ifelse(GalStuff$Highlander_output$best == 'CMA', 1, 2)
      iter = GalStuff$Highlander_output$iteration
      
      ## radio predictions
      flux_3GHz = GalStuff$ProSpectSEDlike_Output$SEDout$FinalFlux$flux[684]
      flux_1.4GHz = GalStuff$ProSpectSEDlike_Output$SEDout$FinalFlux$flux[731]
      
      ## Half age
      magemax = GalStuff$Data$arglist$magemax
      halfmass = GalStuff$ProSpectSEDlike_Output$SEDout$Stars$masstot/2e9
      LookbackTime = cosdistTravelTime(GalStuff$Data$arglist$z, ref = 'Planck15')*1e9
      
      
      massfunc_snorm_trunc_age = function (age, mSFR = 10, mpeak = 10, mperiod = 1, mskew = 0.5, 
                                           magemax = 13.8, mtrunc = 2){
        temp = massfunc_snorm_trunc(age*1e9, mSFR, mpeak, mperiod, mskew, magemax) * age*1e9
        invisible(temp)
      }
      
      values = toppost[1,]
      age = (integral(massfunc_snorm_trunc_age, xmin = 0, xmax = agemax/1e9, 
                      mSFR = 10^values$mSFR,
                      mpeak = values$mpeak,
                      mperiod = 10^values$mperiod,
                      mskew = values$mskew,
                      magemax = (agemax / 1e9))/values$masstot)
      
      
      
      test = c(GalID, Redshift, toppost[1,], BinMasses['TotSMstar'], rangeofvals[,1], TotSMstarlimits[1], rangeofvals[,2], TotSMstarlimits[2] , AGNfrac, type, iter, flux_3GHz, flux_1.4GHz,age)
      test = unlist(test)
      return(test)
   }
  
  Sample = data.table(Sample)
  
  
  parm.names=c('mSFR','mpeak','mperiod','mskew','tau_birth','tau_screen','alpha_SF_birth','alpha_SF_screen','Zfinal','AGNlum', 'AGNan', 'AGNta')
  mon.names=c("LP","masstot","dustmass.birth", "dustmass.screen", "dustmass.total", "dustlum.birth", "dustlum.screen", 
              "dustlum.total", "SFRburst", paste("flux.",filters,sep=''))
  
  names(Sample) = c('CATAID', 'z', parm.names, mon.names, 'StellarMass', paste0(c(parm.names,mon.names,'StellarMass'), '_LB'), 
                    paste0(c(parm.names,mon.names,'StellarMass'), '_UB'), 
                    'AGNfrac', 'Type', 'Iteration', 'flux_3GHz', 'flux_1.4GHz','Age')
  write.csv(Sample, file = paste0(OutputDir,  version_no, '/MasterCat_BestFitParams4.csv'), row.names = F)

}

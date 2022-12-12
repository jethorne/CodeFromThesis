#Loading all the required packages 
#.libPaths(c("/group/pawsey0160/software/sles12sp3/apps/sandybridge/gcc/4.8.5/r/3.6.3/lib64/R/library",.libPaths()))
.libPaths(c("/group/pawsey0160/jthorne/R/magnus/",.libPaths()))
library(ProSpect)
library(magicaxis)
library(data.table)
library(foreach)
library(doParallel)
CoreNumber = 24
# library(sn)
library(LaplacesDemon)
library(celestial)
library(cmaeshpc)
library(scales)
library(Highlander)
registerDoParallel(cores=CoreNumber)

#Setting the fitting version number
version_no = 'v3.42/FIR'

Directory = paste("/group/pawsey0160/jthorne/DEVILS/ProSpect/D10/Outputs/", version_no, '/',  sep='')
OutputDir = Directory
print(paste('output directory:',Directory, sep=''))

if(!dir.exists(OutputDir)){
  dir.create(OutputDir)
}

source('/group/pawsey0160/jthorne/DEVILS/ProSpect/D10/Code/ProSpect_Functions.R')
source('/group/pawsey0160/jthorne/DEVILS/ProSpect/D10/Code/ProSpectSEDlike_plotting.R')

Comb = fread('/group/pawsey0160/jthorne/DEVILS/ProSpect/D10/Data/FIR_Galaxies2fit.csv')
Comb$UID = as.character(format(Comb$UID, scientific = F))


#v3.40 getinfo -- keeping the PAH bands
getinfo=function(CAT=574689, z=0.1, cenwave = cenwave){
  out=Comb[UID==CAT,list(UID, flux_FUV, flux_NUV, flux_u, flux_g,flux_r,flux_i,flux_z, flux_Y, flux_J, flux_H, flux_Ks, flux_S36, flux_S45, flux_S58, flux_S80, flux_M24, flux_M70, flux_P100, flux_P160, flux_S250, flux_S350, flux_S500,
                         flux_err_FUV, flux_err_NUV, flux_err_u,flux_err_g,flux_err_r,flux_err_i,flux_err_z,flux_err_Y, flux_err_J, flux_err_H, flux_err_Ks, flux_err_S36, flux_err_S45, flux_err_S58, flux_err_S80, flux_err_M24, flux_err_M70, flux_err_P100, flux_err_P160,flux_err_S250, flux_err_S350, flux_err_S500)]
  flux=out[,list(flux_FUV, flux_NUV, flux_u, flux_g,flux_r,flux_i,flux_z, flux_Y, flux_J, flux_H, flux_Ks, flux_S36, flux_S45, flux_S58, flux_S80, flux_M24, flux_M70, flux_P100, flux_P160, flux_S250, flux_S350, flux_S500)]
  fluxerr=out[,list( flux_err_FUV, flux_err_NUV, flux_err_u,flux_err_g,flux_err_r,flux_err_i,flux_err_z,flux_err_Y, flux_err_J, flux_err_H, flux_err_Ks, flux_err_S36, flux_err_S45, flux_err_S58, flux_err_S80, flux_err_M24, flux_err_M70, flux_err_P100, flux_err_P160,flux_err_S250, flux_err_S350, flux_err_S500)]
  flux=cbind(flux=as.numeric(flux),fluxerr=as.numeric(fluxerr))
  
  FIRflux = flux[16:22,]
  FUV2NIR = flux[1:15,]
  
  FIRflux[is.na(FIRflux[,'flux']) & !is.na(FIRflux[,'fluxerr']), 'flux'] = 0
 

  flux = rbind(FUV2NIR, FIRflux)
  
  ## removing bands in the PAH features
  #flux[which(cenwave > (1+z)*5e4 & cenwave < (1+z)*15e4),] = NA
  
  floor = rep(0.1, 22)
  for(ii in 1:length(fluxerr)){
    fluxerr[[ii]]=sqrt(fluxerr[[ii]]^2+(flux[,'flux'][[ii]]*floor[ii])^2)
  }
 
  flux=cbind(flux=as.numeric(flux[,'flux']),fluxerr=as.numeric(fluxerr))
 
  flux[flux<0]=0
  data_table = data.table(filters, cenwave)
  flux=cbind(data_table, flux)
  return(list(Lv01=out,flux=flux))
}




#Loading in the required spec librariess and filters for fitting
data(BC03lr)
data(Dale_Msol)
data(EAZY_filters)
filters = c('galex_fuv', 'galex_nuv', 'cfht_u', 'hsc_g','hsc_r','hsc_i','hsc_z', 
            'vista_Y', 'vista_J', 'vista_H', 'vista_Ks','irac_1', 'irac_2', 
            'irac_3', 'irac_4', 'mips_24', 'mips_70', 'pacs_100', 'pacs_160', 
            'spire_250', 'spire_350', 'spire_500')

## Due to weirdness in the names of the VISTA filters it's currently a generic string that works with grep for all filters - this could still end in catastrophe
EAZY_filter_names = c("CAPAK/galex1500.res","CAPAK/galex2500.res", "COSMOS/CFHT_filter_u.txt", "g_HSC.txt","r_HSC.txt"  , 
                      "i_HSC.txt", "z_HSC.txt" , 'VISTA/.s?_system', "IRAC/irac_tr1_2004-08-09.dat", "IRAC/irac_tr2_2004-08-09.dat", 
                      "IRAC/irac_tr3_2004-08-09.dat" ,"IRAC/irac_tr4_2004-08-09.dat", "mips/24" , "mips/70",  "herschel/pacs/100" , "herschel/pacs/160",
                      "herschel/spire/200", "herschel/spire/350", "herschel/spire/500")

EAZY_filter_locations = foreach(ii = 1:length(EAZY_filter_names), .combine = 'c')%do%{
  loc = grep(EAZY_filter_names[ii], EAZY_filters$info, ignore.case = T)
  return(loc)
}

filtout= foreach(jj = 1:length(EAZY_filter_locations))%do%{approxfun(EAZY_filters$filters[EAZY_filter_locations[jj]][[1]])}

cenwave=foreach(i =1:22, .combine = 'c')%do%{cenwavefunc(EAZY_filters$filters[EAZY_filter_locations[i]][[1]])}


# The prior function used
pfunc=function(parm){
  dnorm(parm['tau_birth'],mean=-0.2,sd=0.5,log=TRUE)+
    (-20*erf(parm['tau_screen']-2))+
    dnorm(parm['alpha_SF_birth'],mean=2,sd=1,log=TRUE)+
    dnorm(parm['alpha_SF_screen'],mean=2,sd=1,log=TRUE)+
    (100*erf(parm['mperiod']+2) - 100)
}


IDs = read.csv(commandArgs(trailingOnly=TRUE)[1])
IDs$CATAID = as.character(format(IDs$CATAID ,  scientific = F))


registerDoParallel(cores=CoreNumber)

galaxy_parms = foreach(ii=1:length(IDs$CATAID), .combine='rbind')%dopar%{
  GalID = as.character(format(IDs$CATAID[ii], scientific = F))
  GalID = gsub(" ", "", GalID, fixed = TRUE)
  
  seed = as.numeric(GalID) %% 100 + 1
  
  Redshift = Comb[Comb$UID == GalID, ]$zBest
  print(GalID)
  
  #Setting mpeak limitis 
  LookbackTime = cosdistTravelTime(Redshift, ref = 'Planck15')*1e9
  agemax = 13.38e9 - LookbackTime
  upperlimit = (agemax/1e9) 
  lowerlimit = (0 - LookbackTime)/1e9
  
  Data=list(flux=getinfo(GalID, z = Redshift, cenwave =cenwave )$flux,
            arglist=list(
              z=Redshift,
              massfunc=massfunc_snorm_trunc,
              Z=Zfunc_massmap_lin,
              LumDist_Mpc = cosdistLumDist(z = Redshift, ref = 'Planck15'),  
              agemax=agemax,
              Zagemax = (agemax/1e9) ,
              magemax=(agemax/1e9),
              AGNrm = 60,
              AGNbe = -0.5, 
              AGNal = 4.0,
              AGNct = 100
            ),
            speclib=BC03lr, 
            Dale=Dale_NormTot,
            filtout=filtout,
            SFH=SFHfunc,
            AGN = Fritz, 
            verbose=FALSE,
            Dale_M2L_func=Dale_M2L_func, 
            parm.names=c('mSFR','mpeak','mperiod','mskew','tau_birth','tau_screen','alpha_SF_birth','alpha_SF_screen','Zfinal','AGNlum', 'AGNan', 'AGNta'), 
            mon.names=c("LP","masstot","dustmass.birth", "dustmass.screen", "dustmass.total", "dustlum.birth", "dustlum.screen", 
                        "dustlum.total", "SFRburst", paste("flux.",filters,sep='')), 
            logged=c(T,F,T,F,T,T,F,F,T,T, F, T),
            intervals=list(lo=c(-3,-2, log10(0.3),-0.5,-2.5,-5,0,0,-4, 35, 0.001, -1), 
                           hi=c(4,upperlimit,2,1,1.5,1,4,4,-1.3,49, 89.990, 1)), 
            fit = 'LD', 
            N=length(filters),
            prior=pfunc
  )
  
  Data$flux$cenwave = cenwave
  
  startpoint = (Data$intervals$lo+Data$intervals$hi)/2
  startpoint[10] = 44
  testHigh = Highlander(startpoint, Data,
                        ProSpectSEDlike, Niters=c(2000,2000),  NfinalMCMC = 2000, lower=Data$intervals$lo,
                        upper=Data$intervals$hi, seed=seed, optim_iters = 2, likefunctype = 'LD')

  
  Data$fit = 'check'
  bestfit=ProSpectSEDlike(testHigh$par, Data=Data)
  
  pdf(file = paste0(OutputDir, GalID, '_AGN.pdf' ), width = 11, height = 8)
  plot(bestfit)
  dev.off()
  
  
  ## removing some things to make the outputs smaller
  bestfit$SEDout$StarsAtten = NULL
  bestfit$SEDout$StarsUnAtten = NULL
  #bestfit$SEDout$DustEmit = NULL
  #bestfit$SEDout$Stars$wave_lum = NULL
  #bestfit$SEDout$Stars$lum_unatten = NULL
  #bestfit$SEDout$Stars$lum_atten = NULL
  bestfit$SEDout$Stars$lumtot_unatten = NULL
  bestfit$SEDout$Stars$lumtot_atten = NULL
  bestfit$SEDout$Stars$lumtot_birth = NULL
  bestfit$SEDout$Stars$lumtot_screen = NULL
  
  bestfit$Data$AGN = NULL
  
  bestfit$Data$Dale = NULL
  bestfit$Data$speclib = NULL
  
  testHigh$LD_last$Call = NULL
  testHigh$LD_last$Model = NULL
  
  
  Data$speclib = NULL
  Data$Dale = NULL
  Data$AGN = NULL
  
  print('making output')
  output = list(Data = Data, Highlander_output = testHigh, 
                ProSpectSEDlike_Output = bestfit
  )
  
  print('saving outputs')
  saveRDS(output, file = paste0(OutputDir, GalID, '.rds'))
  print('done!')
  
  }

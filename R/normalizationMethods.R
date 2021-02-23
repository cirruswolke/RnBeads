########################################################################################################################
## normalizationMethods.R
## created: 2021-02-04
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## Implementation of the normalization methods.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

# reimplementation of the basic scaling normalization method for dye bias correction
#
# method='reference', default Illumina scaling normalization to a reference sample, taken from methylumi
# method='internal', singe-sample correlction, taken from minfi
#

rnb.norm.scaling<-function(rnb.set, method=c("internal", "reference")[1L], ref.sample=NULL){
    
    if(!method %in% c("reference", "internal")){
        rnb.error("Unsupported value for method")
    }
    
    qc_int<-qc(rnb.set)
    ctrl_target<-gsub("probes", "controls", rnb.set@target)
    qc_annot<-rnb.get.annotation(ctrl_target, rnb.set@assembly)
    annot<-annotation(rnb.set)
    
    channels<-c("Cy5", "Cy3")
    tokens<-get.platform.tokens(rnb.set@target)
    
    get_norm_int<-function(ch){
        ids<-qc_annot[[tokens$id_col]][grep(tokens$norm_token[ch], qc_annot[[tokens$trg_col]])]
        ids<-as.character(ids)
        qc_int[[ch]][ids,]
    }
    
    norm_int_list<-lapply(channels, get_norm_int)
    names(norm_int_list)<-channels
    
    #taken from methylumi
    Grn.avg <- colMeans(norm_int_list[["Cy3"]])
    Red.avg <- colMeans(norm_int_list[["Cy5"]])
    R.G.ratio = Red.avg/Grn.avg
    if(method=="reference"){
        if(is.null(ref.sample)) ref.sample <- which.min( abs(R.G.ratio-1) )
        rnb.info(paste('Using sample number', ref.sample, 'as reference level...'))
        
        ref <- (Grn.avg + Red.avg)[ref.sample]/2
        if(is.na(ref)) rnb.error("'reference' refers to an array that is not present")
        Grn.factor <- ref/Grn.avg
        Red.factor <- ref/Red.avg
    }else if(method=="internal"){
        Grn.factor <- rep(1, length(R.G.ratio))
        Red.factor <- 1/R.G.ratio
    }
    
    probe.categories<-INTENSITY.SUMMARIZATION.INFO[[rnb.set@target]]

    if(!is.null(rnb.set@status) && rnb.set@status$disk.dump){
        M_temp<-convert.to.ff.matrix.tmp(matrix(NA_real_, nrow=dim(rnb.set@M)[1], ncol=dim(rnb.set@M)[2]))
        U_temp<-convert.to.ff.matrix.tmp(matrix(NA_real_, nrow=dim(rnb.set@U)[1], ncol=dim(rnb.set@U)[2]))
    }else{
        M_temp<-matrix(NA_real_, nrow=nrow(rnb.set@M), ncol=ncol(rnb.set@M))
        U_temp<-matrix(NA_real_, nrow=nrow(rnb.set@U), ncol=ncol(rnb.set@U))
    }
    
    for(pc in probe.categories){
        
        if(pc$Color=="Both"){
            M.factor<-Grn.factor
            U.factor<-Red.factor
        }else if(pc$Color=="Red"){
            M.factor<-Red.factor
            U.factor<-Red.factor
        }else if(pc$Color=="Grn"){
            M.factor<-Grn.factor
            U.factor<-Grn.factor
        }
        indices<-which(annot$Design==pc$Design & annot$Color==pc$Color)
        
        M_temp[indices,]<-sweep(rnb.set@M[indices,], 2, FUN="*", M.factor)
        U_temp[indices,]<-sweep(rnb.set@U[indices,], 2, FUN="*", U.factor)
        
    }
    
    if(!is.null(rnb.set@status) && rnb.set@status$disk.dump && isTRUE(rnb.set@status$discard.ff.matrices)){
        delete(rnb.set@M)
        delete(rnb.set@U)
    }
    
    rnb.set@M<-M_temp
    rnb.set@U<-U_temp
    
    if(!is.null(rnb.set@status) && rnb.set@status$disk.dump){
        rm(M_temp); rm(U_temp); rnb.cleanMem()
    }
    
    rnb.set<-update.meth(rnb.set)
    
    return(rnb.set)
}



# reimplementation of the basic subtraction  normalization method

rnb.bgcorr.subtr<-function(rnb.set, bg.quant=0.05, offset=0){
    
    qc_int<-qc(rnb.set)
    ctrl_target<-gsub("probes", "controls", rnb.set@target)
    qc_annot<-rnb.get.annotation(ctrl_target, rnb.set@assembly)
    annot<-annotation(rnb.set)
    
    channels<-c("Cy5", "Cy3")
    tokens<-get.platform.tokens(rnb.set@target)
    
    get_bg_int<-function(ch){
        ids<-qc_annot[[tokens$id_col]][grep(tokens$bg_token, qc_annot[[tokens$trg_col]])]
        ids<-as.character(ids)
        qc_int[[ch]][ids,]
    }
    
    bg_int_list<-lapply(channels, get_bg_int)
    names(bg_int_list)<-channels
    
    bg<-list()
    for(ch in channels) {
        bg[[ch]]<-colQuantiles(bg_int_list[[ch]], probs=bg.quant)
    }
    
    probe.categories<-INTENSITY.SUMMARIZATION.INFO[[rnb.set@target]]
    
    if(!is.null(rnb.set@status) && rnb.set@status$disk.dump){
        M_temp<-convert.to.ff.matrix.tmp(matrix(NA_real_, nrow=dim(rnb.set@M)[1], ncol=dim(rnb.set@M)[2]))
        U_temp<-convert.to.ff.matrix.tmp(matrix(NA_real_, nrow=dim(rnb.set@U)[1], ncol=dim(rnb.set@U)[2]))
    }else{
        M_temp<-matrix(NA_real_, nrow=nrow(rnb.set@M), ncol=ncol(rnb.set@M))
        U_temp<-matrix(NA_real_, nrow=nrow(rnb.set@U), ncol=ncol(rnb.set@U))
    }
    
    for(pc in probe.categories){
        
        if(pc$Color=="Both"){
            bg.M<-bg[["Cy5"]]
            bg.U<-bg[["Cy3"]]
        }else if(pc$Color=="Red"){
            bg.M<-bg[["Cy3"]]
            bg.U<-bg[["Cy3"]]
        }else if(pc$Color=="Grn"){
            bg.M<-bg[["Cy5"]]
            bg.U<-bg[["Cy5"]]
        }
        indices<-which(annot$Design==pc$Design & annot$Color==pc$Color)
        for(i in 1:ncol(rnb.set@M)) {
            M_temp[indices,i] <- pmax(rnb.set@M[indices,i] - bg.M[i], 1)
            U_temp[indices,i] <- pmax(rnb.set@U[indices,i] - bg.U[i], 1)
        }
    }
    
    if(!is.null(rnb.set@status) && rnb.set@status$disk.dump && isTRUE(rnb.set@status$discard.ff.matrices)){
        delete(rnb.set@M)
        delete(rnb.set@U)
    }
    
    rnb.set@M<-M_temp
    rnb.set@U<-U_temp
    
    if(!is.null(rnb.set@status) && rnb.set@status$disk.dump){
        rm(M_temp); rm(U_temp); rnb.cleanMem()
    }
    
    rnb.set<-update.meth(rnb.set)
    
    return(rnb.set)
    
}


get.platform.tokens<-function(platform){
    
    channels<-c("Cy5", "Cy3")
    token<-setNames(rep(NA_character_, 2), channels)
    
    dict<-list()
    
    if(rnb.set@target=="probes27"){
        dict$norm_token["Cy5"] <- 'Norm.G'
        dict$norm_token["Cy3"] <- 'Norm.R'
        dict$bg_token<-"Negative"
        dict$id_col<-"ID"
        dict$trg_col<-"Target"
    } else if(rnb.set@target=="probes450"){
        dict$norm_token["Cy5"] <- 'Norm_(C|G)'
        dict$norm_token["Cy3"] <- 'Norm_(A|T)'
        dict$bg_token<-"Negative"
        dict$id_col<-"ID"
        dict$trg_col<-"Target"
    } else if(rnb.set@target=="probesEPIC"){
        dict$norm_token["Cy5"] <- 'NORM_(C|G)'
        dict$norm_token["Cy3"] <- 'NORM_(A|T)'
        dict$bg_token<-"NEGATIVE"
        dict$id_col<-"ID"
        dict$trg_col<-"Target"
    }else if(rnb.set@target=="probesMMBC"){
        dict$norm_token["Cy5"] <- 'NORM_(C|G)'
        dict$norm_token["Cy3"] <- 'NORM_(A|T)'
        dict$bg_token<-"NEGATIVE"
        dict$id_col<-"ID"
        dict$trg_col<-"Target"
    }
    
    return(dict)
}


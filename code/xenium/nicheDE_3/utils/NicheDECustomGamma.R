# gamma is cell-type specific
niche_DE_no_parallel = function(object,C = 150,M = 10,print = T,Int = T,batch = T,self_EN = F, gamma){
  print("Running custom NicheDE inference")  
  nb_lik = function(x,mu,disp){
    #returns negative log likelihood: Var = mu + mu^2/size
    return(-sum(dnbinom(x=x, size = disp, mu = mu,log = TRUE)))
  }
  niche_DE_core = function(counts,constants,iter){
    #get reference expression vector
    ref_expression = constants$ref_expr[,iter]
    #get expected expression of gene j
    EEJ = as.vector(constants$num_cells%*%ref_expression)
    #get gene number
    if(iter%%1000 == 0){
      print(paste0("Gene ",iter," out of ",constants$gene_total))
    }

    valid = 0
    liks_val = NA
    n_type = ncol(constants$effective_niche)
    #check if we need to do niche-DE
    if((sum(counts)<constants$C) | (mean(ref_expression<constants$CT_filter)!=1)==F){
      null = c(1:n_type^2)
      liks_val = NA
    }
    if((sum(counts)>constants$C)&(mean(ref_expression<constants$CT_filter)!=1)){
      #get pstg matrix
      pstg = constants$num_cells%*%as.matrix(diag(ref_expression))/EEJ
      pstg[,ref_expression<constants$CT_filter] = 0
      pstg[pstg<0.05]=0
      #get X
      X = matrix(NA,nrow(pstg),n_type^2)
      #get counter
      for(k in c(1:nrow(pstg))){
        #get spatial covariates
        X[k,] = as.vector(round(constants$effective_niche[k,],2)%*%t(as.matrix(pstg[k,])))
        if(constants$self_EN == F){
          ind_remove = (c(1:n_type)-1)*n_type + c(1:n_type)
          X[k,ind_remove] = 0
        }
      }
      rm("pstg")
      #gc()
      #get index, niche pairs that are non existent
      null = which(apply(X,2,function(x){sum(x>0,na.rm = T)})<constants$M)
      X_partial = X
      rest = c(1:ncol(X))
      if(length(null)>0){
        X_partial = X[,-null]
        rest = rest[-null]
      }
      rm("X")
      #add batch
      if(constants$batch == T &length(unique(constants$batchID)) > 1){
        batchvar = as.factor(constants$batchID)
      }else{
        batchvar = rep(1,length(counts))
      }
      #gc()
      #continue if at least one index,niche pair is viable
      if(length(null)!=n_type^2 & constants$Int == T){
        tryCatch({
          nvar = ncol(X_partial)
          #if expected expression for a spot is 0, remove it
          bad_ind  = which(EEJ==0)
          #run neg binom regression
          if(length(bad_ind)>0){
            full_glm =suppressWarnings({glm(counts[-bad_ind]~X_partial[-bad_ind,] + batchvar[-bad_ind]+ offset(log(EEJ[-bad_ind])), family = "poisson")}) #do full glm
          }else{
            full_glm = suppressWarnings({glm(counts~X_partial + batchvar + offset(log(EEJ)), family = "poisson")}) #do full glm
          }
          mu_hat = as.numeric(exp(predict(full_glm)))#get mean
          #get dicpersion parameter
          A = optimize(nb_lik,x = counts,mu = mu_hat, lower = 0.05, upper = 100) #get overdispersion parameter
          #save dispersion parameter
          disp = A$minimum
          #save likelihood
          liks_val = -A$objective
          #calculate W matrix for distribution of beta hat
          W =as.vector(mu_hat/(1 + mu_hat/disp))#get W matrix
          rm("mu_hat")
          #perform cholesky decomp for finding inverse of X^TWX
          if(length(bad_ind)>0){
            X_partial = X_partial[-bad_ind,]
            #X_partial = Matrix::Matrix(X_partial, sparse=TRUE)
          }
          if(constants$batch == T &length(unique(constants$batchID)) > 1){
            if(length(bad_ind)>0){
              dummy_col = fastDummies::dummy_cols(batchvar[-bad_ind],remove_first_dummy = T, remove_selected_columns = T)
            }else{
              dummy_col = fastDummies::dummy_cols(batchvar,remove_first_dummy = T, remove_selected_columns = T)
            }
            #append dummy variables and intercept
            X_partial = cbind(X_partial,dummy_col)
            X_partial = cbind(X_partial,rep(1,nrow(X_partial)))
          }else{
            #append dummy variables and intercept
            X_partial = cbind(X_partial,rep(1,nrow(X_partial)))
          }
          X_partial = as.matrix(X_partial)
          X_partial = Matrix::Matrix(X_partial, sparse=TRUE)
          #get variance matrix. Variance is [t(X_partial*W)%*%X_partial]^(-1)
          var_mat = Matrix::t(X_partial*W)%*%X_partial
          rm("X_partial")
          #gc()
          #if there are degenerate columns, remove them
          new_null = c()

          #account for batch ID variables
          new_null = which(diag(as.matrix(var_mat))==0)
          if(length(new_null)>0){
            var_mat = var_mat[-new_null,-new_null]
            null = sort(c(null,rest[new_null[new_null <= nvar]]))
            new_null = new_nul[new_null <= var]
          }

          if(length(null)!=n_type^2){
            #get coefficients for important variables and intercept
            coeff = full_glm$coefficients[1:(nvar+1)]
            coeff[is.na(coeff)] = 0
            #cholesky decomposition
            A = Matrix::chol(var_mat,LDL = FALSE,perm = FALSE)
            #get covaraince matrix
            V = Matrix::solve(A)%*%Matrix::t(Matrix::solve(A))
            V = V[1:(n_type^2-length(null)),1:(n_type^2-length(null))]
            #save V as an upper triangular matrix
            V =  Matrix::Matrix(upper.tri(V,diag = T)*V, sparse=TRUE)
            #remove large objects
            rm("A")
            rm("var_mat")
            #get sd matrix
            tau = sqrt(Matrix::diag(V))
            V_ = matrix(NA,n_type,n_type)
            if(length(null)==0){
              V_ = matrix(tau,n_type,n_type)
            }else{
              V_[c(1:n_type^2)[-null]] = tau}
            #get beta
            beta = matrix(NA,n_type,n_type)
            if(length(new_null)>0){
              beta[c(1:n_type^2)[-null]] = coeff[-c(1,new_null+1)]
            }

            if(length(new_null)==0){
              if(length(null)==0){
                beta = matrix(coeff[-c(1)],n_type,n_type)
              }else{
                beta[c(1:n_type^2)[-null]] = coeff[-c(1)]}
            }
            #record test statitistic
            T_ = beta/V_

            beta = Matrix::t(beta)
            T_ = Matrix::t(T_)
            #make rownames and column names correspond to cell types
            colnames(beta) = paste0("n: ",constants$cell_types)
            rownames(beta) = paste0("i: ",constants$cell_types)
            colnames(T_) = paste0("n: ",constants$cell_types)
            rownames(T_) = paste0("i: ",constants$cell_types)
            #record that niche-de worked
            valid = 1
          }
          #end of if statement
        }, #get pval
        error = function(e) {
          #skip_to_next <<- TRUE
          })
      }
      if(length(null)!=n_type^2 & constants$Int == F){
        tryCatch({
          nvar = ncol(X_partial)
          #if expected expression for a spot is 0, remove it
          lm = lm((counts - EEJ) ~ X_partial + batchvar)
          rm("X_partial")
          #gc()
          sum_lm =  summary(lm)
          #get log likelihood
          liks_val = stats::logLik(lm)
          #get vcov mat
          V = (sum_lm$cov.unscaled)*(sum_lm$sigma^2)
          #see if any betas have 0 variance
          new_null = which(diag(as.matrix(V))==0)
          if(length(new_null)>0){
            V = V[-new_null,-new_null]
            null = sort(c(null,rest[new_null[new_null <= nvar]]))
            new_null = new_null[new_null <= nvar]
          }
          #remove first observation
          V = V[-c(1),-c(1)]
          #save V as an upper triangular matrix
          V =  Matrix::Matrix(upper.tri(V,diag = T)*V, sparse=TRUE)
          #get beta coefficients
          if(length(null)!=n_type^2){
            #get coefficients
            coeff = lm$coefficients[1:(nvar+1)]
            coeff[is.na(coeff)] = 0
            #get number of true coefficients
            num_coef = n_type^2 - length(null)
            coeff_T = sum_lm$coefficients[c(1:(num_coef+1)),3]
            coeff_T[is.na(coeff_T)] = 0
            #get beta vector
            beta = matrix(NA,n_type,n_type)
            T_ = matrix(NA,n_type,n_type)
            if(length(new_null)>0){
              beta[c(1:n_type^2)[-null]] = coeff[-c(1,new_null+1)]
              T_[c(1:n_type^2)[-null]] = coeff_T[-c(1)]
            }

            if(length(new_null)==0){
              if(length(null)==0){
                beta = matrix(coeff[-c(1)],n_type,n_type)
                T_ = matrix(coeff_T[-c(1)],n_type,n_type)
              }else{
                beta[c(1:n_type^2)[-null]] = coeff[-c(1)]
                T_[c(1:n_type^2)[-null]] = coeff_T[-c(1)]
              }
            }
            #record test statistic
            beta = Matrix::t(beta)
            T_ = Matrix::t(T_)
            #make rownames and column names correspond to cell types
            colnames(beta) = paste0("n: ",constants$cell_types)
            rownames(beta) = paste0("i: ",constants$cell_types)
            colnames(T_) = paste0("n: ",constants$cell_types)
            rownames(T_) = paste0("i: ",constants$cell_types)
            #record that niche-DE worked
            valid = 1
          }
          #end of if statement
        }, #get pval
        error = function(e) {
          skip_to_next <<- TRUE})
      }
    }

    if(iter == constants$gene_total){
      print("Niche-DE Complete")
    }
    if(valid == 1){
      A = list(T_stat = T_,betas = beta, Varcov = V,nulls = null, valid = valid,log_likelihood = liks_val)
      rm(list=ls()[! ls() %in% c("A")])
      #gc()
      return (A)
    }else{
      A = list(T_stat = NA,betas = NA, Varcov = NA,nulls = c(1:n_type^2), valid = 0,log_likelihood = liks_val)
      rm(list=ls()[! ls() %in% c("A")])
      #gc()
      return (A)
    }

  }

  #starting Message
  print(paste0('Starting Niche-DE analysis with parameters C = ',C,', M = ',M,', gamma = '))
  print(gamma)
  #initialize list output
  object@niche_DE = vector(mode = "list", length = length(object@sigma))
  names(object@niche_DE) = object@sigma
  #get counter: which bandwidth we are on
  counter = 1
  #get validity matrix
  valid = matrix(0,ncol(object@counts),length(object@sigma))
  #iterate over kernel bandwidths
  for(sig in object@sigma){
    print(paste0('Performing Niche-DE analysis with kernel bandwidth:',sig,' (number ',counter,' out of ',length(object@sigma),' values)'))
    #get expression filter (gamma)
    CT_filter = sapply(rownames(object@ref_expr), function(ct){
        quantile(object@ref_expr[ct, ], gamma[ct])
    })
    names(CT_filter) = rownames(object@ref_expr)
    #CT_filter = apply(object@ref_expr,1,function(x){quantile(x,gamma)})
    #get number of genes
    ngene = ncol(object@counts)
    #save nuber of cell types
    n_type = ncol(object@num_cells)
    #save dimnames for arrays
    dimnames = list(A = colnames(object@num_cells),B  = colnames(object@num_cells), C = colnames(object@counts))
    dimdims = c(ncol(object@num_cells),ncol(object@num_cells),ncol(object@counts))

    #initialize constant parameter list
    constant_param = list(effective_niche = object@effective_niche[[counter]], num_cells = object@num_cells,
                          CT_filter = CT_filter,C = C,M = M,Int = Int,gene_total = ngene,batch = batch,
                          batchID = object@batch_ID,ref_expr = object@ref_expr,cell_types = object@cell_types,self_EN = self_EN)

    print(paste0("Starting Niche-DE"))
    results = vector("list", length = ngene)

    for(j in c(1:ngene)){
      results[[j]] = niche_DE_core(counts = object@counts[,j],constants = constant_param,iter = j)
    }

    names(results) = object@gene_names
    #save valid matrix
    valid[,counter] = unlist(lapply(results, function(result) result$valid))
    #save object
    object@niche_DE[[counter]] = results
    counter = counter + 1
    #remove everything but what's needed
    print("Cleaning disk for next iteration")
    rm(list=ls()[! ls() %in% c("object","counter","nb_lik","niche_DE_core","C","M","gamma","valid","cores","Int","batch","print","nicheDE","self_EN")])
    gc()
  }
  #get column sums of counts matrix to see how many genes pass filtering
  A = rowSums(valid)
  #get number of genes that pass filtering
  num_pass = sum(A>=1,na.rm = T)
  print('Computing Positive Niche-DE Pvalues')
  object = get_niche_DE_pval_fisher(object,pos = T)
  print('Computing Negative Niche-DE Pvalues')
  object = get_niche_DE_pval_fisher(object,pos = F)
  print(paste0('Niche-DE analysis complete. Number of Genes with niche-DE T-stat equal to ',num_pass))
  if(num_pass < 1000){
    warning('Less than 1000 genes pass. This could be due to insufficient read depth of data or size of C parameter. Consider changing choice of C parameter')
  }
  return(object)

}

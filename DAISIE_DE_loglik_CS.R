
DAISIE_DE_loglik_CS <- function(
    pars1 = pars1,
    pars2 = pars2,
    datalist = datalist,
    methode = methode,
    rtol,
    atol
)
{
  cond = pars2[3]
  if(length(pars1) == 5)
  {
    logp0 = DAISIE_DE_logp0 (datalist, pars1, methode)
    if(is.null(datalist[[1]]$not_present))
    {
      loglik = (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) * logp0
      numimm = (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) + length(datalist) - 1
    } else {
      loglik = datalist[[1]]$not_present * logp0
      numimm = datalist[[1]]$not_present + length(datalist) - 1
    }
    logcond = (cond == 1) * log(1 - exp(numimm * logp0))
    for(i in 2:length(datalist))
    {
      datalist[[i]]$type1or2 = 1
    }
  } else {
    numimm = datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2 + length(datalist) - 1
    numimm_type2 = length(which(unlist(datalist)[which(names(unlist(datalist)) == "type1or2")] == 2))
    numimm_type1 = length(datalist) - 1 - numimm_type2
    if(is.na(pars1[11]) == FALSE)
    {
      if(pars1[11] < numimm_type2/numimm | pars1[11] > (1 - numimm_type1 /numimm)) { return(-Inf) }
      datalist[[1]]$not_present_type2 = max(0,round(pars1[11] * numimm) - numimm_type2)
      datalist[[1]]$not_present_type1 = numimm - (length(datalist) - 1) - datalist[[1]]$not_present_type2
    }
    logp0_type1 = DAISIE_DE_logp0 (datalist, pars1, methode)
    logp0_type2 = DAISIE_DE_logp0 (datalist, pars1, methode)
    loglik = datalist[[1]]$not_present_type1 * logp0_type1 + datalist[[1]]$not_present_type2 * logp0_type2
    logcond = (cond == 1) * log(1 - exp((datalist[[1]]$not_present_type1 + numimm_type1) * logp0_type1 + (datalist[[1]]$not_present_type2 + numimm_type2) * logp0_type2))
  }
  loglik = loglik - logcond
  vec_likelihood <- c()
  for (i in 2:length(datalist))
  {
    if (datalist[[i]]$stac == 1 )
    {
      likelihood <- DAISIE_DE_logpNE_max_age_coltime(datalist, i, pars1, methode,rtol, atol)
      vec_likelihood <- c(vec_likelihood, likelihood)
      s1 = sprintf('Status of colonist: %d, Parameters: %f %f %f %f %f',1,pars1[1],pars1[2],pars1[3],pars1[4],pars1[5])
      s2 = sprintf(', Loglikelihood: %f',likelihood)
      cat(s1,s2,"\n",sep = "")
      
    }
    else if (datalist[[i]]$stac == 2 && length(datalist[[i]]$branching_times) == 2 && datalist[[i]]$missing_species == 0)
      
    {
      likelihood <-  DAISIE_DE_logpES(datalist, i, pars1, methode,rtol, atol)
      vec_likelihood <- c(vec_likelihood, likelihood)
      s3 = sprintf('Status of colonist: %d, Parameters: %f %f %f %f %f',2,pars1[1],pars1[2],pars1[3],pars1[4],pars1[5])
      s4 = sprintf(', Loglikelihood: %f',likelihood)
      cat(s3,s4,"\n",sep = "")
    }
    else if (datalist[[i]]$stac == 3 && length(datalist[[i]]$branching_times) == 2 && datalist[[i]]$missing_species == 0)
      
    {
      likelihood <-  DAISIE_DE_logpES_mainland(datalist, i, pars1, methode,rtol, atol) 
      vec_likelihood <- c(vec_likelihood, likelihood)
      s3 = sprintf('Status of colonist: %d, Parameters: %f %f %f %f %f',2,pars1[1],pars1[2],pars1[3],pars1[4],pars1[5])
      s4 = sprintf(', Loglikelihood: %f',likelihood)
      cat(s3,s4,"\n",sep = "")
    }
    else if (datalist[[i]]$stac == 3 && length(datalist[[i]]$branching_times) > 2)
      
    {
      likelihood <-  DAISIE_DE_logpEC_mainland(datalist, i, pars1, methode, rtol, atol ) 
      vec_likelihood <- c(vec_likelihood, likelihood)
      s3 = sprintf('Status of colonist: %d, Parameters: %f %f %f %f %f',2,pars1[1],pars1[2],pars1[3],pars1[4],pars1[5])
      s4 = sprintf(', Loglikelihood: %f',likelihood)
      cat(s3,s4,"\n",sep = "")
    }
    else if (datalist[[i]]$stac == 4 )
    {
      likelihood <- DAISIE_DE_logpNE(datalist, i, pars1, methode, rtol, atol)
      vec_likelihood <- c(vec_likelihood, likelihood)
      s5 = sprintf('Status of colonist: %d, Parameters: %f %f %f %f %f', 4, pars1[1],pars1[2],pars1[3],pars1[4],pars1[5])
      s6 = sprintf(', Loglikelihood: %f',likelihood)
      cat(s5,s6,"\n",sep = "")
    }
    else if (datalist[[i]]$stac == 5 )
    {
      likelihood <- DAISIE_DE_logpES_max_age_coltime(datalist, i, pars1, methode, rtol, atol)
      vec_likelihood <- c(vec_likelihood, likelihood)
      s7 = sprintf('Status of colonist: %d, Parameters: %f %f %f %f %f', 5 ,pars1[1],pars1[2],pars1[3],pars1[4],pars1[5])
      s8 = sprintf(', Loglikelihood: %f',likelihood)
      cat(s7,s8,"\n",sep = "")
    }
    else if (datalist[[i]]$stac == 6 && length(datalist[[i]]$branching_times) > 2 )
    {
      likelihood <- function_Factor_Loglik_EC_max_Age_approximation1(datalist, i, pars1, methode,rtol, atol)
      vec_likelihood <- c(vec_likelihood, likelihood)
      s9 = sprintf('Status of colonist: %d, Parameters: %f %f %f %f %f', 6, pars1[1], pars1[2],pars1[3],pars1[4],pars1[5])
      s10 = sprintf(', Loglikelihood: %f',likelihood)
      cat(s9,s10,"\n",sep = "")
    }
    else if (datalist[[i]]$stac == 2 && length(datalist[[i]]$branching_times) > 2)
    {
      
      likelihood <- DAISIE_DE_logpEC(datalist, i, pars1, methode, rtol, atol) 
      vec_likelihood <- c(vec_likelihood, likelihood)
      s11 = sprintf('Status of colonist: %d, Parameters: %f %f %f %f %f',2,pars1[1],pars1[2],pars1[3],pars1[4],pars1[5])
      s12 = sprintf(', Loglikelihood: %f',likelihood)
      cat(s11,s12,"\n",sep = "")
    }
    else if (datalist[[i]]$stac == 2 && length(datalist[[i]]$branching_times) == 2 && datalist[[i]]$missing_species != 0)
    {
      
      likelihood <- function_Factor_Loglik_EC_approximation20(datalist, i, pars1) + function_DDD_Likelihood_EC(data_list, i, pars1)
      vec_likelihood <- c(vec_likelihood, likelihood)
      s11 = sprintf('Status of colonist: %d, Parameters: %f %f %f %f %f',2,pars1[1],pars1[2],pars1[3],pars1[4],pars1[5])
      s12 = sprintf(', Loglikelihood: %f',likelihood)
      cat(s11,s12,"\n",sep = "")
    }
    else if (datalist[[i]]$stac == 6 && length(datalist[[i]]$branching_times) == 2 && datalist[[i]]$missing_species != 0)
    {
      
      likelihood <- function_Factor_Loglik_EC_max_Age_approximation_2(datalist, i, pars1)
      vec_likelihood <- c(vec_likelihood, likelihood)
      s11 = sprintf('Status of colonist: %d, Parameters: %f %f %f %f %f',6,pars1[1],pars1[2],pars1[3],pars1[4],pars1[5])
      s12 = sprintf(', Loglikelihood: %f',likelihood)
      cat(s11,s12,"\n",sep = "")
    }
    #new_likelihood <- sum(vec_likelihood)
    
  }
  ss = sprintf(' Parameter s: %f %f %f %f %f',pars1[1],pars1[2],pars1[3],pars1[4],pars1[5])
  ss2 = sprintf(', Loglikelihood: %f',sum(vec_likelihood)+ loglik)
  cat(ss,ss2,"\n",sep = "")
  return (sum(vec_likelihood) + loglik)
}

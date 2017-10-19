
runninginCluster = 0

if(runninginCluster == 1){
  args <- commandArgs(TRUE)
  dataset = args[1]
  id = args[2]
  inFile = paste0('/nfs/guille/bugid/adams/emmott_backup/bake_off/data/mothersets/binary/', dataset, '/', dataset, '.preproc.csv')
  RDSFilePath = '/nfs/cluster-fserv/cluster-share/amran/RRFAnalyst/particle/'
  explDir = paste0('/nfs/cluster-fserv/cluster-share/amran/all.6.explanations/', dataset, '/explanations/')
  refidxDir = paste0('/nfs/cluster-fserv/cluster-share/amran/refIdx/', dataset, '/')
  analystRespFolder = paste0('/nfs/cluster-fserv/cluster-share/amran/all.6.explanations/', dataset, '/ResponseRRFAnalyst/')
  HRnkDir <- paste0('/nfs/cluster-fserv/cluster-share/amran/all.6.explanations/', dataset, '/')
}else{
  # set current directory
  directory = "."
  dataset = 'shuttle'
  id = 100
  inFile = paste0(directory, '/', dataset, '/', dataset, '.preproc.csv')
  RDSFilePath = paste0(directory, '/', dataset, '/')
  explDir = paste0(directory, '/', dataset, '/explanations/')
  refidxDir = paste0(directory, '/', dataset, '/refidx/')
  analystRespFolder = paste0(directory, '/', dataset,'/ResponseRRFAnalyst/')
  HRnkDir <- paste0(directory, '/', dataset, '/')
}

# ------------------------------------------------------------------

allProbs <- readRDS(paste0(RDSFilePath, 'allModel.Probs.', dataset, '.RDS'))

X <- read.csv(inFile, header = TRUE)
anIndx <- which(X$ground.truth=='anomaly')
y <- X$ground.truth
X$ground.truth = NULL

MAX_FEAT <- 6
getAnalystResp <- function(idx, rnk){
  val <- rep(0, ncol(X))
  mask = 0
  for(i in 1:MAX_FEAT){
    val[i] <- 1
    mask = mask + 2^(ncol(X)-rnk[i])
    tryCatch({val[i] <- allProbs[idx, paste0(mask)]},
             error=function(cond) {print(c(mask, 'model does not exist'))})
  }
  return (val)
}


refFile = paste0('refidx_', id, '.csv')
# read corresponding anomaly reference indices in ground truth explanation
refIdx = read.csv(paste0(refidxDir, refFile))

"HRnkIdx <- read.csv(paste0(HRnkDir, 'HRankBenchIdx.csv'))
HRnkIdx <- HRnkIdx[which(HRnkIdx$benchId == id), ]
benchAnoId <- rep(0, nrow(HRnkIdx))
motherAnoId <- rep(0, nrow(HRnkIdx))
c <- 1
for(hidx in HRnkIdx$anoIdx){
  t <- which(refIdx$anoIdx == hidx)
  benchAnoId[c] <- t
  motherAnoId[c] <- which(anIndx == refIdx$refIdx[t])
  c <- c + 1
}"

benchAnoId <- rep(0, nrow(refIdx))
motherAnoId <- rep(0, nrow(refIdx))

c <- 1
for(anom in refIdx$refIdx)
{
  benchAnoId[c] <- c
  motherAnoId[c] <- which(anIndx == anom)
  c <- c + 1
}

# print(motherAnoId)


explanation.seq_dropout   = read.csv(paste0(explDir, 'exp_', dataset, '_', id, '_seq_dropout.csv'))
explanation.seq_marginal  = read.csv(paste0(explDir, 'exp_', dataset, '_', id, '_seq_marginal.csv'))

maxRow <- 1000
result.seq_dropout <- matrix(0, nrow = maxRow, ncol = ncol(X)+3)
result.seq_marginal <- matrix(0, nrow = maxRow, ncol = ncol(X)+3)
cnt <- 1

numAnoInst = length(motherAnoId)
for(ii in 1:numAnoInst){
  rank.seq_dropout  = as.numeric(explanation.seq_dropout[ii, -1])
  rank.seq_marginal = as.numeric(explanation.seq_marginal[ii, -1])

  i <- benchAnoId[ii]
  b_id <- refIdx$anoIdx[i]
  r_id <- refIdx$refIdx[i]
  m_id <- motherAnoId[ii]
  result.seq_dropout[cnt, ]  <- c(id, b_id, r_id, getAnalystResp(m_id, rank.seq_dropout))
  result.seq_marginal[cnt, ] <- c(id, b_id, r_id, getAnalystResp(m_id, rank.seq_marginal))
  cnt <- cnt + 1
}

writeToFile <- function(result, cnt, fName){
  colnames(result) <- c('benchID', 'anoIdx', 'refIdx', paste0('Resp.upto.', 1:ncol(X)))
  write.csv(as.data.frame(result)[1:(cnt-1), ], paste0(analystRespFolder,fName), row.names = F, quote = F)
}

analyst = 'RRF'
writeToFile(result = result.seq_dropout,  cnt, paste0(analyst, 'Analyst.Resp.benchExpl.', id, '.seq_dropout.', dataset, '.csv'))
writeToFile(result = result.seq_marginal, cnt, paste0(analyst, 'Analyst.Resp.benchExpl.', id, '.seq_marginal.', dataset, '.csv'))

# -----------------------------------------------------------------

expMethods = c('seq_dropout', 'seq_marginal')

getMFP <- function(resp, TH){
  consumed = rep(0, nrow(resp))
  for(i in 1:nrow(resp)){
    cnt = 1
    for(j in 1:ncol(resp)){
      if(resp[i, j] <= TH)
        break;
      cnt = cnt + 1
    }
    consumed[i] = cnt
  }
  return (consumed)
}


respFilePrefix = paste0(directory, '/', dataset, '/ResponseRRFAnalyst/RRFAnalyst.Resp.benchExpl.')

BenchIDST = 100
BenchIDEnd = 100
if(dataset == 'concrete' || dataset == 'skin')
  BenchIDEnd = 1200

response = list()
cnt = 0
assigned = rep(0, 4)
idx = rep(0, 1000)
for(id in BenchIDST:BenchIDEnd){
  if(id %% 10 == 0)
    print(id)
  for(i in 1:4){
    fName = paste0(respFilePrefix, id, '.', expMethods[i], '.', dataset, '.csv')
    if(file.exists(fName)){
      tempResp = read.csv(fName)
      if(assigned[i] == 0){
        response[[i]] = tempResp
        assigned[i] = 1
      }else{
        response[[i]] = rbind(response[[i]], tempResp)
      }
    }
  }
}
  
result <- matrix(0, nrow = 3*nrow(response[[1]]), ncol = 4)
colnames(result) <- c('benchID', 'anoIdx', 'refIdx', 'MFP')
for(i in 1:2){
  result[,1] <- rep(response[[1]]$benchID, 3)
  result[,2] <- rep(response[[1]]$anoIdx, 3)
  result[,3] <- rep(response[[1]]$refIdx, 3)
  result[,4] <- c(getMFP(resp = response[[i]][,-(1:3)], TH = 0.1),
                  getMFP(resp = response[[i]][,-(1:3)], TH = 0.2),
                  getMFP(resp = response[[i]][,-(1:3)], TH = 0.3))
  write.csv(result, paste0('MFP.',expMethods[i], '.', dataset, '.allBenchR10.csv'), row.names = F, quote = F)
}

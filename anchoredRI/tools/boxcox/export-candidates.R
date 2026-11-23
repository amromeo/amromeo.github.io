# Summaries are derived from the saved per-dataset candidate tables, before any safeguard.
bc_export_tables <- function(directory, include_arms=FALSE) {
  t <- RIbench::loadTestsetDefinition(); t <- t[t$N<=50000,]
  grid <- expand.grid(replicate=1:20,case=1:45)
  meta <- rbind(data.frame(key=paste0('r-',t$Index),stage=ifelse(t$N<=5000,'small','50000'),
                           analyte=t$Analyte),
                data.frame(key=sprintf('s-%02d-%02d',grid$case,grid$replicate),stage='simulation',
                           analyte=stress_design$family[grid$case]))
  stopifnot(nrow(meta)==5220L,!anyDuplicated(meta$key))
  write_table <- function(frame,file) write.csv(frame,file.path(directory,file),row.names=FALSE)
  for (method in c('anchoredRI',if(include_arms)c('armA','armB'))) {
    checkpoint <- if(method=='anchoredRI')'validation' else 'fits'
    pairs <- lapply(meta$key,function(key) {
      r <- readRDS(file.path(directory,checkpoint,paste0(key,'.rds')))
      if(method!='anchoredRI')r <- r[[method]]
      stopifnot(all(r$candidate_table$key==key),nrow(r$candidate_table)>0L)
      list(candidates=r$candidate_table,gap=r$bic_gap)
    })
    candidates <- do.call(rbind,lapply(pairs,'[[','candidates'))
    gaps <- cbind(meta,do.call(rbind,lapply(pairs,'[[','gap'))[,-1,drop=FALSE])
    write_table(candidates,paste0(method,'-candidates.csv'))
    write_table(gaps,paste0(method,'-bic-gaps.csv'))
    if(method=='anchoredRI') {
      write_table(gaps[which(gaps$within_2),],'baseline-near-tie-limits.csv')
      groups <- c(list(all=seq_len(nrow(gaps))),split(seq_len(nrow(gaps)),gaps$stage),
                  split(seq_len(nrow(gaps)),paste(gaps$stage,gaps$analyte,sep='/')))
      summary <- do.call(rbind,lapply(names(groups),function(group)
        cbind(group=group,bc_near_tie_counts(gaps[groups[[group]],,drop=FALSE]))))
      write_table(summary,'baseline-near-tie-summary.csv')
    }
  }
}

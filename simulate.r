library(dplyr)
library(simulateGP)
library(glue)
library(data.table)
library(ggplot2)
library(tidyr)

bfile_dir <- "~/repo/opengwas-api-internal/opengwas-api/app/ld_files"
pops <- c("AFR", "AMR", "EAS", "EUR", "SAS")

for(pop in pops){
    glue("plink2 --bfile {bfile_dir}/{pop} --freq --out {pop}") %>% system()
}

bim <- lapply(pops, \(pop) fread(file.path(bfile_dir, paste0(pop, ".bim"))) %>% mutate(pop = pop))
names(bim) <- pops

freqs <- lapply(pops, \(pop) fread(paste0(pop, ".afreq")) %>% mutate(pos = bim[[pop]]$V4, pop=pop))
names(freqs) <- pops


# Choose 10k SNPs to make the polygenic component
# Choose 1 SNP to have a large effect that only has an effect in individuals who drink alcohol.
# GWAS on 10k SNPs to identify EUR discovery variants


o <- lapply(freqs, \(x) {
    x %>% filter(ALT_FREQS > 0.2 & ALT_FREQS < 0.8) %>% select(ID, ALT_FREQS, pop, REF, ALT)
}) %>% bind_rows()

common_snps <- o %>% group_by(ID) %>% summarise(n=n()) %>% filter(n == 5)

set.seed(12345)

rsids <- sample(common_snps$ID, 10000, replace=FALSE)
map <- freqs[[1]] %>% select(snp=ID, chr="#CHROM", pos=pos, ea=ALT, oa=REF, af=ALT_FREQS) %>% filter(snp %in% rsids)

intrsid <- rsids[1]

params <- map %>% generate_gwas_params(h2=0.5, S=-0.4, Pi=1)
params$beta[params$snp == intrsid] <- 0.035

params2 <- params
params2$beta[params2$snp == intrsid] <- 0.035 * 0.1
temp <- subset(freqs[["SAS"]], ID %in% map$snp) %>% select(snp=ID, af=ALT_FREQS)
params2 <- params2 %>% select(-c(af))
params2 <- inner_join(params2, temp, by="snp")
params


ss <- lapply(pops, \(pop) {
    if(pop == "SAS"){
        ss <- generate_gwas_ss(params2, nid=ifelse(pop == "EUR", 100000, 30000)) %>% mutate(pop=pop)
    } else {
        p <- params %>% select(-c(af))
        temp <- subset(freqs[[pop]], ID %in% map$snp) %>% select(snp=ID, af=ALT_FREQS)
        p <- inner_join(p, temp, by="snp")
        ss <- generate_gwas_ss(p, nid=ifelse(pop == "EUR", 100000, 30000)) %>% mutate(pop=pop)
    }
    ss
})
names(ss) <- pops

lapply(ss, \(x) x %>% filter(pval < 5e-8) %>% nrow())

eur_disc <- ss[[4]] %>% filter(pval < 5e-8)
eur_disc %>% filter(snp == intrsid)


head(params)

write.table(params %>% select(snp, ea, beta), "score.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(params2 %>% select(snp, ea, beta), "score2.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

scores <- lapply(pops, \(pop) {
    if(pop == "SAS") {
        glue("plink2 --bfile {bfile_dir}/{pop} --score score2.txt --out {pop}") %>% system()
    } else {
        glue("plink2 --bfile {bfile_dir}/{pop} --score score.txt --out {pop}") %>% system()
    }
    score <- fread(paste0(pop, ".sscore"))
    score %>% rename(score=SCORE1_AVG) %>% mutate(score=score * 1000, pop=pop)
})

scores <- bind_rows(scores)
group_by(scores, pop) %>% summarise(m=mean(score), s=sd(score))

ggplot(scores, aes(x=score, fill=pop)) + geom_density(alpha=0.5)
ggplot(scores, aes(x=pop, y=score)) + geom_violin()


write.table(params %>% filter(snp %in% eur_disc$snp) %>% select(snp, ea, beta), "score.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(params2 %>% filter(snp %in% eur_disc$snp) %>% select(snp, ea, beta), "score2.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)


scores2 <- lapply(pops, \(pop) {
    # write.table(ss[[pop]] %>% filter(snp %in% eur_disc$snp) %>% select(snp, ea, bhat), "score.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
    if(pop == "SAS") {
        glue("plink2 --bfile {bfile_dir}/{pop} --score score2.txt --out {pop}") %>% system()
    } else {
        glue("plink2 --bfile {bfile_dir}/{pop} --score score.txt --out {pop}") %>% system()
    }
    score <- fread(paste0(pop, ".sscore"))
    score %>% rename(score=SCORE1_AVG) %>% mutate(score=score * 1000, pop=pop, method="eur_disc")
}) %>% bind_rows()


group_by(scores2, pop) %>% summarise(m=mean(score), s=sd(score))
ggplot(scores2, aes(x=score, fill=pop)) + geom_density(alpha=0.5)
ggplot(scores2, aes(x=pop, y=score)) + geom_violin() + geom_boxplot(width=0.2)

## Heterogeneity of effect estimates

gwas_res <- lapply(pops, \(pop) {
    sel <- ss[[pop]] %>% filter(snp %in% eur_disc$snp) %>%
        select(snp, ea, oa, af, bhat, se, pval)
    names(sel)[4:7] <- paste0(names(sel)[4:7], "_", pop)
    sel
}) %>%
    Reduce(inner_join, .)


a <- tidyr::pivot_longer(gwas_res %>% select(snp, starts_with("af_")), cols=c("af_EUR", "af_AFR", "af_EAS", "af_AMR", "af_SAS"), names_to="pop", values_to="af") %>% mutate(pop=gsub("af_", "", pop)) 

merge(a, a, by="snp") %>% ggplot(aes(x=af.x, y=af.y, colour=pop.x)) + geom_point() + geom_abline() + geom_density_2d(alpha=0.5) + facet_wrap(pop.x~pop.y)


table(gwas_res$pval_EUR < 5e-8)
table(gwas_res$pval_AFR < 5e-8)
table(gwas_res$pval_EAS < 5e-8)
table(gwas_res$pval_SAS < 5e-8)
table(gwas_res$pval_AMR < 5e-8)

table(gwas_res$pval_EUR < 0.005)
table(gwas_res$pval_AFR < 0.005)
table(gwas_res$pval_EAS < 0.005)
table(gwas_res$pval_SAS < 0.005)
table(gwas_res$pval_AMR < 0.005)




fixed_effects_meta_analysis <- function(beta_vec, se_vec) {
    w <- 1 / se_vec^2
    beta <- sum(beta_vec * w, na.rm=T) / sum(w, na.rm=T)
    se <- sqrt(1 / sum(w, na.rm=T))
    pval <- pnorm(abs(beta / se), lower.tail = FALSE)
    Qj <- w * (beta-beta_vec)^2
    Q <- sum(Qj, na.rm=T)
    Qdf <- sum(!is.na(beta_vec))-1
    if(Qdf == 0) Q <- 0
    Qjpval <- pchisq(Qj, 1, lower.tail=FALSE)
    Qpval <- pchisq(Q, Qdf, lower.tail=FALSE)
    return(list(beta=beta, se=se, pval=pval, Q=Q, Qdf=Qdf, Qpval=Qpval, Qj=Qj, Qjpval=Qjpval))
}

for(i in 1:nrow(gwas_res)) {
    a <- fixed_effects_meta_analysis(gwas_res[i, c("bhat_EUR", "bhat_AFR", "bhat_EAS", "bhat_AMR", "bhat_SAS")] %>% unlist(),
                                gwas_res[i, c("se_EUR", "se_AFR", "se_EAS", "se_AMR", "se_SAS")] %>% unlist())
    gwas_res$Qpval[i] <- a$Qpval
}


table(gwas_res$Qpval < 0.05/nrow(gwas_res))




het <- gwas_res[gwas_res$Qpval < 0.05/nrow(gwas_res), ]

# Make long format

inner_join(
    tidyr::pivot_longer(het %>% select(snp, starts_with("bhat")), cols=c("bhat_EUR", "bhat_AFR", "bhat_EAS", "bhat_AMR", "bhat_SAS"), names_to="pop", values_to="bhat") %>% mutate(pop=gsub("bhat_", "", pop)),
    tidyr::pivot_longer(het %>% select(snp, starts_with("se")), cols=c("se_EUR", "se_AFR", "se_EAS", "se_AMR", "se_SAS"), names_to="pop", values_to="se") %>% mutate(pop=gsub("se_", "", pop))
) %>% 
    ggplot(aes(x=bhat, y=pop, colour=pop)) + 
        geom_point() +
        geom_errorbarh(aes(xmin=bhat-1.96*se, xmax=bhat+1.96*se), height=0) +
        geom_vline(xintercept=0, linetype="dashed")



# GxE interaction


g <- lapply(pops, \(pop) {
    glue("plink2 --bfile {bfile_dir}/{pop} --snps {intrsid} --recode A --out {pop}") %>% system()
    g <- fread(paste0(pop, ".raw"))
    g$pop <- pop
    g
}) %>% bind_rows()

dat <- tibble(
    id=g$FID,
    true_score = scores$score,
    eur_score = scores2$score,
    rs58903867 = g[,7] %>% unlist(),
    pop=g$pop
)

dat$alcohol <- sample(c(0, 1), nrow(dat), replace=TRUE, prob=c(0.2, 0.8))
dat$alcohol[dat$pop == "SAS"] <- sample(c(0, 1), sum(dat$pop == "SAS"), replace=TRUE, prob=c(0.9, 0.1))
dat$phen <- dat$true_score + dat$rs58903867 * dat$alcohol * 0.3 + rnorm(nrow(dat), 0, 1)
dat




summary(lm(phen ~ rs58903867, data=dat))

summary(lm(phen ~ rs58903867 * alcohol, data=dat))
summary(lm(phen ~ rs58903867 * pop, data=dat))

group_by(dat, pop) %>% do({
    a <- summary(lm(phen ~ rs58903867, data=.))$coefficients
 tibble(bhat=a[2,1], se=a[2,2], pval=a[2,4])
}) %>%
    ggplot(aes(x=bhat, y=pop, colour=pop)) + 
        geom_point() +
        geom_errorbarh(aes(xmin=bhat-1.96*se, xmax=bhat+1.96*se), height=0) +
        geom_vline(xintercept=0, linetype="dashed")


group_by(dat, alcohol) %>% do({
    a <- summary(lm(phen ~ rs58903867, data=.))$coefficients
 tibble(bhat=a[2,1], se=a[2,2], pval=a[2,4])
}) %>%
    ggplot(aes(x=bhat, y=as.factor(alcohol), colour=alcohol)) + 
        geom_point() +
        geom_errorbarh(aes(xmin=bhat-1.96*se, xmax=bhat+1.96*se), height=0) +
        geom_vline(xintercept=0, linetype="dashed")


write.csv(gwas_res, file="gwas_res.csv")
write.csv(dat, file="individual_values.csv")








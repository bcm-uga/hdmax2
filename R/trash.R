# Charger le package nécessaire
library(MASS)

# Définir les paramètres
set.seed(123)  # Pour la reproductibilité
n <- 100  # Nombre d'observations
p = 200
p1 <- p/4   # Nombre de variables dans le premier paquet
p2 <- p/4   # Nombre de variables dans le deuxième paquet
p3 <- p/4  
p4 <- p/4 # Nombre de variables dans le troisième paquet
total_vars = p   # Total de variables

# Fonction pour générer un bloc avec une corrélation fixe pour toutes les paires
generate_fixed_corr_block <- function(size, corr_range = c(0.5, 0.8)) {
  corr_val <- runif(1, min = corr_range[1], max = corr_range[2])  # Générer une corrélation fixe
  block <- matrix(corr_val, nrow = size, ncol = size)  # Remplir la matrice avec cette corrélation
  diag(block) <- 1  # Mettre 1 sur la diagonale
  return(block)
}

# Générer des matrices de corrélation pour chaque paquet
corr_block1 <- generate_fixed_corr_block(p1)
corr_block2 <- generate_fixed_corr_block(p2)
corr_block3 <- generate_fixed_corr_block(p3)
corr_block4 <- generate_fixed_corr_block(p4)

# Construire la matrice de corrélation complète avec corrélation zéro entre les paquets
corr_matrix <- matrix(0, nrow = total_vars, ncol = total_vars)
corr_matrix[1:p1, 1:p1] <- corr_block1
corr_matrix[(p1+1):(p1+p2), (p1+1):(p1+p2)] <- corr_block2
corr_matrix[(p1+p2+1):total_vars, (p1+p2+1):total_vars] <- corr_block3
corr_matrix[(p1+p2+p3+1):total_vars, (p1+p2+p3+1):total_vars] <- corr_block4
# Vérifier si la matrice est semi-définie positive (nécessaire pour mvrnorm)
is_positive_definite <- function(mat) {
  return(all(eigen(mat)$values > 0))
}

if (!is_positive_definite(corr_matrix)) {
  stop("La matrice de corrélation n'est pas positive définie. Essayez de réajuster.")
}

# Générer des données avec la corrélation spécifiée
data <- mvrnorm(n = n, mu = rep(0, total_vars), Sigma = corr_matrix)

# Convertir les données en data frame pour une meilleure lisibilité
colnames(data) <- paste0("Var", 1:total_vars)
data <- as.data.frame(data)

# Afficher un aperçu des données générées
head(data)





prop.causal.probes = 0.1
n=200
p=500

## effects/correlation
rho = 0.5


sigma = 0.2
sd.A = 0.2
mean.A = 0.8
mean.B = 0.8
mean.A = 0.6
sd.B = 0.2
mean.B = 0.8













library(timereg)

data(sTRACE)
# Fits Aalen model 
out<-aalen(Surv(time,status==9)~age+sex+diabetes+chf+vf,
           sTRACE,max.time=7,n.sim=100)

summary(out)
par(mfrow=c(2,3))
plot(out)

# Fits semi-parametric additive hazards model 
out<-aalen(Surv(time,status==9)~const(age)+const(sex)+const(diabetes)+chf+vf,
           sTRACE,max.time=7,n.sim=100)

summary(out)
par(mfrow=c(2,3))
plot(out)

## Excess risk additive modelling 
data(mela.pop)
dummy<-rnorm(nrow(mela.pop));

# Fits Aalen model  with offsets 
out<-aalen(Surv(start,stop,status==1)~age+sex+const(dummy),
           mela.pop,max.time=7,n.sim=100,offsets=mela.pop$rate,id=mela.pop$id,
           gamma=0)
summary(out)
par(mfrow=c(2,3))
plot(out,main="Additive excess riks model")

# Fits semi-parametric additive hazards model  with offsets 
out<-aalen(Surv(start,stop,status==1)~age+const(sex),
           mela.pop,max.time=7,n.sim=100,offsets=mela.pop$rate,id=mela.pop$id)
summary(out)
plot(out,main="Additive excess riks model")


library(survival)
cancer = data(cancer)

mod_aalen_leukemia = timereg::aalen(Surv(leukemia$time,leukemia$status)~leukemia$x)
mod_aalen_leukemia$pval.testBeq0

mod_aalen_lung = timereg::aalen(Surv(lung$time,lung$status)~lung$ph.ecog)
mod_aalen_lung$obs.testBeq0

mod_aalen_colon = timereg::aalen(Surv(colon$time,colon$status)~colon$nodes+ colon$age+colon$sex)
mod_aalen_colon$pval.testBeq0[2]




)

mod_aalen_colon$table

sum = summary(out)$Test
sum[1,]
summary(out)[1,]
out$pval.testBeq0[2]

coef_cum <- out$cum[nrow(out$cum), -1]  # Dernière ligne, sans la colonne de temps

# Extraire la matrice de variance robuste des coefficients (robvar.cum)
var_robust <- out$robvar.cum[, nrow(out$robvar.cum)]  # On prend la dernière période

# Calcul des erreurs standards
se_robust <- sqrt(diag(var_robust))

# Calculer les statistiques de Wald
wald_stats <- coef_cum / se_robust

# Calcul des p-values associées (bilatérales)
p_values_precises <- 2 * pnorm(-abs(wald_stats))